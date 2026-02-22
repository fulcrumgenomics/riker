use proc_macro::TokenStream;
use quote::{ToTokens, format_ident, quote};
use syn::{Data, DeriveInput, Fields, Lit, Meta, parse_macro_input};

/// Derive macro that extracts `///` doc comments from a struct and its fields,
/// generating an implementation of `crate::metrics::MetricDocs`.
///
/// Field names default to the Rust identifier but respect `#[serde(rename = "...")]`.
///
/// # Panics
/// Panics if applied to a non-struct item or a struct without named fields.
#[proc_macro_derive(MetricDocs)]
pub fn derive_metric_docs(input: TokenStream) -> TokenStream {
    let input = parse_macro_input!(input as DeriveInput);
    let name = &input.ident;

    // Extract struct-level doc comments.
    let struct_doc = extract_doc_string(&input.attrs);

    // Extract per-field docs.
    let fields = match &input.data {
        Data::Struct(data) => match &data.fields {
            Fields::Named(named) => &named.named,
            _ => panic!("MetricDocs only supports structs with named fields"),
        },
        _ => panic!("MetricDocs can only be derived for structs"),
    };

    let mut field_names = Vec::new();
    let mut field_descs = Vec::new();

    for field in fields {
        let field_name = serde_rename(field)
            .unwrap_or_else(|| field.ident.as_ref().expect("named field").to_string());
        let field_doc = extract_doc_string(&field.attrs);
        field_names.push(field_name);
        field_descs.push(field_doc);
    }

    let n = field_names.len();

    let expanded = quote! {
        impl crate::metrics::MetricDocs for #name {
            fn metric_description() -> &'static str {
                #struct_doc
            }

            fn field_docs() -> &'static [crate::metrics::FieldDoc] {
                static DOCS: [crate::metrics::FieldDoc; #n] = [
                    #(
                        crate::metrics::FieldDoc {
                            name: #field_names,
                            description: #field_descs,
                        },
                    )*
                ];
                &DOCS
            }
        }
    };

    TokenStream::from(expanded)
}

/// Join all `#[doc = "..."]` attributes into a single trimmed string.
fn extract_doc_string(attrs: &[syn::Attribute]) -> String {
    let lines: Vec<String> = attrs
        .iter()
        .filter_map(|attr| {
            if !attr.path().is_ident("doc") {
                return None;
            }
            if let Meta::NameValue(nv) = &attr.meta {
                if let syn::Expr::Lit(expr_lit) = &nv.value {
                    if let Lit::Str(lit) = &expr_lit.lit {
                        return Some(lit.value());
                    }
                }
            }
            None
        })
        .collect();

    lines.iter().map(|l| l.trim()).collect::<Vec<_>>().join(" ").trim().to_string()
}

/// Attribute macro that passes through the original Options struct unchanged and
/// directly generates a `Multi{OriginalName}` struct with `--{prefix}::field-name`
/// CLI args, plus `validate()` and `From<OriginalName>` conversions.
///
/// Fields are classified into three categories:
/// - **Optional**: `Option<T>` fields are kept as `Option<T>` in the Multi struct.
/// - **Defaulted**: Non-`Option<T>` fields with `default_value_t` keep their type and default.
/// - **Required**: Non-`Option<T>` fields without `default_value_t` are wrapped in `Option<T>`
///   in the Multi struct. The generated `validate()` method checks they are `Some`.
///
/// Usage: `#[multi_options("prefix", "Help Heading")]`
///
/// # Panics
/// Panics if applied to a non-struct, a struct without named fields, or a struct
/// containing `#[command(flatten)]` (structs must be flat).
#[proc_macro_attribute]
pub fn multi_options(attr: TokenStream, item: TokenStream) -> TokenStream {
    // Parse attribute arguments: ("prefix", "heading")
    let args = syn::parse::<MultiOptionsArgs>(attr).expect(
        "multi_options requires two string literal arguments: #[multi_options(\"prefix\", \"heading\")]"
    );
    let prefix = &args.prefix;
    let heading = &args.heading;

    let input = syn::parse::<syn::ItemStruct>(item).expect("multi_options requires a struct");
    let struct_name = &input.ident;
    let multi_struct_name = format_ident!("Multi{}", struct_name);

    // Extract named fields.
    let fields = match &input.fields {
        syn::Fields::Named(named) => &named.named,
        _ => panic!("multi_options only supports structs with named fields"),
    };

    // Validate: no #[command(flatten)] fields.
    for field in fields {
        for attr in &field.attrs {
            if attr.path().is_ident("command") {
                let tokens = attr.meta.to_token_stream().to_string();
                assert!(
                    !tokens.contains("flatten"),
                    "multi_options does not support #[command(flatten)] on field `{}`. \
                     Inline the flattened struct's fields directly.",
                    field.ident.as_ref().unwrap()
                );
            }
        }
    }

    // Build per-field data.
    let mut multi_field_tokens = Vec::new();
    let mut validate_field_tokens = Vec::new();
    let mut from_original_field_tokens = Vec::new();

    for field in fields {
        generate_field_tokens(
            field,
            struct_name,
            prefix,
            &mut multi_field_tokens,
            &mut validate_field_tokens,
            &mut from_original_field_tokens,
        );
    }

    // Emit the original struct unchanged + the Multi struct + validate + From<Original>.
    let original = input.to_token_stream();

    let expanded = quote! {
        #original

        /// Prefixed options struct generated by `#[multi_options]` for the multi command.
        #[derive(clap::Args, Debug, Clone)]
        #[command(next_help_heading = #heading)]
        pub struct #multi_struct_name {
            #(#multi_field_tokens)*
        }

        impl #multi_struct_name {
            /// Validate required fields and convert to the original options struct.
            pub fn validate(self) -> anyhow::Result<#struct_name> {
                let opts = self;
                Ok(#struct_name {
                    #(#validate_field_tokens)*
                })
            }
        }

        impl From<#struct_name> for #multi_struct_name {
            fn from(opts: #struct_name) -> Self {
                Self {
                    #(#from_original_field_tokens)*
                }
            }
        }
    };

    TokenStream::from(expanded)
}

/// Generate Multi-struct field, validate assignment, and From<Original> assignment
/// for a single field of the Options struct.
fn generate_field_tokens(
    field: &syn::Field,
    struct_name: &syn::Ident,
    prefix: &str,
    multi_field_tokens: &mut Vec<proc_macro2::TokenStream>,
    validate_field_tokens: &mut Vec<proc_macro2::TokenStream>,
    from_original_field_tokens: &mut Vec<proc_macro2::TokenStream>,
) {
    let field_ident = field.ident.as_ref().expect("named field");
    let field_type = &field.ty;
    let kebab_name = field_ident.to_string().replace('_', "-");
    let doc_comment = extract_doc_string(&field.attrs);
    let is_option = is_option_type(field_type);
    let has_default = has_default_value_t(field);
    let preserved = extract_preserved_arg_attrs(field);

    let long_name = format!("{prefix}::{kebab_name}");
    let prefixed_ident = format_ident!("{}_{}", prefix.replace('-', "_"), field_ident);

    // A field is "required" if it's not Option<T> and has no default_value_t.
    let is_required = !is_option && !has_default;

    if is_required {
        let required_doc = format!("{doc_comment} Required when {prefix} is selected.");
        let error_msg = format!("--{long_name} is required when {prefix} is selected");

        multi_field_tokens.push(quote! {
            #[doc = #required_doc]
            #[arg(long = #long_name #(#preserved)*)]
            pub #prefixed_ident: Option<#field_type>,
        });
        validate_field_tokens.push(quote! {
            #field_ident: opts.#prefixed_ident
                .ok_or_else(|| anyhow::anyhow!(#error_msg))?,
        });
        from_original_field_tokens.push(quote! {
            #prefixed_ident: Some(opts.#field_ident),
        });
    } else if is_option {
        multi_field_tokens.push(quote! {
            #[doc = #doc_comment]
            #[arg(long = #long_name #(#preserved)*)]
            pub #prefixed_ident: #field_type,
        });
        validate_field_tokens.push(quote! {
            #field_ident: opts.#prefixed_ident,
        });
        from_original_field_tokens.push(quote! {
            #prefixed_ident: opts.#field_ident,
        });
    } else {
        let default_attr = quote! { , default_value_t = #struct_name::default().#field_ident };
        multi_field_tokens.push(quote! {
            #[doc = #doc_comment]
            #[arg(long = #long_name #default_attr #(#preserved)*)]
            pub #prefixed_ident: #field_type,
        });
        validate_field_tokens.push(quote! {
            #field_ident: opts.#prefixed_ident,
        });
        from_original_field_tokens.push(quote! {
            #prefixed_ident: opts.#field_ident,
        });
    }
}

/// Parsed arguments for the `#[multi_options("prefix", "heading")]` attribute.
struct MultiOptionsArgs {
    prefix: String,
    heading: String,
}

impl syn::parse::Parse for MultiOptionsArgs {
    fn parse(input: syn::parse::ParseStream) -> syn::Result<Self> {
        let prefix_lit: syn::LitStr = input.parse()?;
        input.parse::<syn::Token![,]>()?;
        let heading_lit: syn::LitStr = input.parse()?;
        Ok(Self { prefix: prefix_lit.value(), heading: heading_lit.value() })
    }
}

/// Check if a field has `default_value_t` in its `#[arg(...)]` attribute.
fn has_default_value_t(field: &syn::Field) -> bool {
    for attr in &field.attrs {
        if !attr.path().is_ident("arg") {
            continue;
        }
        if let Ok(nested) = attr
            .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        {
            for meta in &nested {
                match meta {
                    Meta::NameValue(nv) if nv.path.is_ident("default_value_t") => return true,
                    Meta::Path(path) if path.is_ident("default_value_t") => return true,
                    _ => {}
                }
            }
        }
    }
    false
}

/// Check if a type is `Option<T>`.
fn is_option_type(ty: &syn::Type) -> bool {
    if let syn::Type::Path(type_path) = ty {
        if let Some(segment) = type_path.path.segments.last() {
            return segment.ident == "Option";
        }
    }
    false
}

/// Extract #[arg(...)] attributes that should be preserved (everything except
/// `long`, `short`, and `default_value_t`), returning them as `, key = value` tokens.
fn extract_preserved_arg_attrs(field: &syn::Field) -> Vec<proc_macro2::TokenStream> {
    let mut preserved = Vec::new();
    for attr in &field.attrs {
        if !attr.path().is_ident("arg") {
            continue;
        }
        if let Ok(nested) = attr
            .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        {
            for meta in &nested {
                match meta {
                    Meta::NameValue(nv) => {
                        let name = nv.path.get_ident().map(ToString::to_string);
                        if !matches!(name.as_deref(), Some("long" | "short" | "default_value_t")) {
                            preserved.push(quote! { , #nv });
                        }
                    }
                    Meta::Path(path) => {
                        let name = path.get_ident().map(ToString::to_string);
                        if !matches!(name.as_deref(), Some("long" | "short")) {
                            preserved.push(quote! { , #path });
                        }
                    }
                    Meta::List(_) => {
                        preserved.push(quote! { , #meta });
                    }
                }
            }
        }
    }
    preserved
}

/// Check for `#[serde(rename = "...")]` and return the rename value if present.
fn serde_rename(field: &syn::Field) -> Option<String> {
    for attr in &field.attrs {
        if !attr.path().is_ident("serde") {
            continue;
        }
        if let Ok(nested) = attr
            .parse_args_with(syn::punctuated::Punctuated::<Meta, syn::Token![,]>::parse_terminated)
        {
            for meta in &nested {
                if let Meta::NameValue(nv) = meta {
                    if nv.path.is_ident("rename") {
                        if let syn::Expr::Lit(expr_lit) = &nv.value {
                            if let Lit::Str(lit) = &expr_lit.lit {
                                return Some(lit.value());
                            }
                        }
                    }
                }
            }
        }
    }
    None
}
