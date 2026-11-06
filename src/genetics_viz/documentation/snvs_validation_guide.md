## SNV Validation Guide

This guide helps you validate genetic variants accurately and consistently.

### 1. Validate the Correct Variant

**Always verify you are validating the variant shown at the top of the dialog.**

- The variant key (chr:pos:ref:alt) is displayed prominently
- The sample ID is shown next to the variant
- Double-check before saving your validation

### 2. Use Precise Inheritance Descriptors

Choose the most specific inheritance pattern that applies:

| Inheritance | When to use |
|-------------|-------------|
| **de novo** | Variant is confirmed absent in both parents |
| **paternal** | Variant inherited from father only |
| **maternal** | Variant inherited from mother only |
| **not paternal** | Variant is absent in father (mother status unknown/unavailable) |
| **not maternal** | Variant is absent in mother (father status unknown/unavailable) |
| **either** | Variant present in both parents |
| **homozygous** | Variant is homozygous in the sample |
| **unknown** | Inheritance cannot be determined |

*e.g.* if we have only the mother and that she carries the variant, we do not put maternal but unknown

### 3. Validation Status

| Status | Meaning |
|--------|---------|
| **present** | Variant is confirmed present in the sample |
| **absent** | Variant is NOT present (false positive call) |
| **uncertain** | Cannot determine if variant is real |
| **different** | A variant exists but differs from what was called |
| **in phase MNV** | Variant is part of a multi-nucleotide variant (confirmed present) |

*different* should be a last resort option and you should write a comment (most likely for indels in LCR region that tries to depict a STR of different length)

### 4. Resolving Conflicts

When multiple validations exist for the same variant:

- Use the **Ignore switch** in the "Previous validations" section
- Toggle ignored validations to exclude them from analysis
- Only non-ignored validations are considered for final status
- Use this to resolve disagreements between validators
- Usually, when 2 users disagree, a third one should resolve the conflict (and add his own validation)

### 5. Comments

Use the comment field to:

- Explain unusual observations
- Note technical issues with the data
- Reference related variants (e.g., MNV components)
- Document reasons for ignoring a validation
- "LCR" is a good comment for repeat regions. On IGV, the reference genome is not capitalized in repeats.

## Tips

- **Add parents** to the IGV viewer to assess inheritance
- **Check coverage** - low coverage regions may have unreliable calls
- **Look at strand bias** - variants on only one strand may be artifacts
- **Consider population frequency** - rare variants need more scrutiny
- **Validate only the variant you are looking at** - not another variant in the same position
