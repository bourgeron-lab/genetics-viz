# Changelog

All notable changes to genetics-viz will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.26.0] - 2026-10-23

### Added
- **Gene and variant pages, across every cohort.** Type a gene symbol or a variant (`chr22:50674000:G:A`, `22-50674000-G-A`) in the new header box, on any page:
  - `/gene/<symbol>` lists the gene's small variants (one row per variant with its carriers, families, cohorts and validation tally, then one row per carrier), WisecondorX calls, STR genotypes past normal, ROH segments and diagnostics, from every cohort of every data directory. A previous symbol (`KIAA1109`) opens the current gene (`BLTP1`) and searches both.
  - `/variant/<variant>` shows who carries it, and every validation and diagnostic ever recorded for it, ignored ones included, even when it is in no cohort table.
  - A family of another data directory opens there: the link switches your data directory first.
- **Gene index.** The pages read a local index of the cohort tables, built at startup and checked every 30 min; only tables that changed are read again. It lives in `gene_index_dir` (default `.gene_index/` next to the config file). Families with no cohort-level table are not in it.

## [0.25.0] - 2026-10-22

### Security
- **Data files need a login.** The data directories (CRAMs, coverage, every TSV) were served to anyone who could reach the server. They now answer 401 without a logged-in session, 404 once a directory is removed on the admin page, and are sent with `Cache-Control: private`. Range requests still work for IGV.
- **Passwords are stored as salted scrypt** instead of unsalted SHA-512. Existing passwords keep working: each user's first login rehashes theirs and saves the config. No reset is needed.
- **Login backs off after repeated failures**: 5 free failures per username and 20 per client address, then a doubling wait of up to 15 min. Unknown usernames are treated like real ones, so the login does not reveal which accounts exist.
- **Validations, diagnostics and wave reviews record the logged-in user** as author, not the value of the readonly User box, which the browser could change.

### Added
- **Shareable links.** The address bar now keeps the view, so it can be sent to a colleague who opens exactly the same view:
  - family pages: the open tab, the selected members, the shown Wombat config with its filters, sort and page, and the SV call types;
  - cohort page: the Workflow tab, the right-hand panel tab and the individuals table's column filters;
  - Search: the whole form, run automatically when the link is opened. Sources are named by display name, so a link carries no server path.

  Only values that differ from the defaults are written, so links stay short (`?tab=roh`). Anything the page cannot show falls back to the default.
- **Login returns to the linked page.** Opening a link while logged out goes to the login, then to that page instead of the home page.

### Changed
- **STR tab: "Hide loci where every member is normal" is on by default.** The tab's Download follows the switch.
- **Rows validated "in phase MNV" count as present** in every validation filter. Every status selection used to hide them, including the default.
- **An empty Impacts selection filters nothing on Search**, as in the Wombat tab. It used to show no rows.
- **All validations** shows a present call's inheritance beyond de novo and homozygous (e.g. `maternal`), as the other pages do.
- A Wombat filter change goes back to the table's first page. Choosing an impact level refreshes the table once instead of once per checkbox.

### Fixed
- **Two curators saving at the same time could corrupt a curation file.** Validations, diagnostics, notes, wave reviews, Sample QC acknowledgements, UPD and STR validations and report conclusions now all go through one locked write. Before:
  - the first two saves to a new file could both write the header;
  - lines reached the file only after the lock was released, so the next writer could find it empty;
  - an ignore toggle or a note deletion emptied the file before locking it, so a save in between could be lost.
- **A note or comment starting with a quote lost it.** A comment containing `"` that went through an ignore toggle or a note deletion before this release was saved quoted, and now shows its quotes doubled.
- **Hidden directories under `cohorts/`** (a NAS `.snapshot`) were listed as incomplete cohorts.
- The family member boxes were always ticked, whatever the selection.

### Refactored
- Validation statuses and their consensus are defined once (`utils/validation_status.py`), and Wombat badges and filters are shared by the Wombat tab, Search and the per-file validation page (`utils/variant_rows.py`).
- The three largest files are split: `sv_dialog.py` (into `sv_details.py`, `sv_igv_views.py`, `sv_validation_panel.py`), `search.py` (`search_form.py`, `search_results.py`, `utils/search_rows.py`) and `wombat_tab.py` (`wombat_rows.py`, `wombat_filters.py`, `wombat_stats.py`).

### Note
- **Behind a proxy, the per-address login limit counts the proxy's address** unless uvicorn trusts the forwarded headers (`FORWARDED_ALLOW_IPS`, with `X-Forwarded-For` passed along). Otherwise every user shares one address, and 20 failures from anyone slow logins for everyone.
- Links keep filter values (and on Search, loci and gene names) in the URL query, so they appear in the server's access log, as family IDs in the paths already do.

## [0.24.1] - 2026-10-21

### Fixed
- **Family report: the findings table ran off the page.** An SNV's consequence was VEP's raw comma list (`downstream_gene_variant,missense_variant,…`), one unbreakable word. It now reads as words, most severe first, without up/downstream and other modifier terms when a coding one is there (`missense variant`). Every report table cell can also wrap now, so no long value widens a table past the page.

## [0.24.0] - 2026-10-21

### Added
- **Family report for the referring clinician.** Both family pages have two new buttons at the end of the title row, for curators and administrators only:
  - **See report** opens the report in a full-screen dialog;
  - **PDF** opens the browser's print dialog. The default file name is `<cohort>_<family>_report_<date>`, and each A4 page carries a confidentiality footer and page numbers.

  The report recaps the family in this order: the pedigree; the curator's conclusion; "Results at a glance"; ancestry, with PCA plots; Sample QC issues and their reviews; runs of homozygosity; four polygenic scores (autism, ADHD, epilepsy, intelligence), each a small distribution with every member's percentile; aneuploidy; UPD; and the curators' pathogenic and uncertain findings, with their validation and inheritance. It ends with short methods and the tool versions, and the pipeline version and commit of the run.
  - Aneuploidy combines the Sample QC karyotype for the sex chromosomes with WisecondorX gains or losses covering most of an autosome.
  - It is addressed to the clinician, not to patients, and says so on every page with a research-use disclaimer.
- **Curator conclusion per family**, written in the report dialog and printed at the top of the report. Saved in `reports/family_summaries.tsv`, one row per save; the latest one counts.

### Changed
- On the cohort family page, the previous/next family arrows now sit at the right end of the title row, after the report buttons.
- `jinja2` is now a direct dependency (it was already installed through NiceGUI).

## [0.23.0] - 2026-10-21

### Added
- **UPD validation and diagnostic.** A uniparental disomy is reviewed from the family ROH tab (a new "UPD candidates" strip under the summary) or from the Sample QC Issues table (the review action on a `UPD_CANDIDATE` flag). Both open the same dialog:
  - the ROH and Sample QC evidence for the child's chromosome;
  - the validation history, and for curators a form: present / absent / uncertain, parental origin (pre-filled from the flag), comment;
  - the diagnostic, through the shared diagnostic dialog.

  Records are keyed by family, child and chromosome (`chr2:UPD`), in `validations/upds.tsv` and `diagnostics/upds.tsv`. UPD is a validation and diagnostic kind, so the All and Statistics pages list it.
- **A validated UPD acknowledges its Sample QC flag**: present as confirmed, absent as explained; uncertain leaves it counting. The header badge, member icons, cohort column and panel, home card and exports follow.

### Changed
- **Family Diagnostics panel: STR and UPD lines say what was found.** An STR reads as its gene badge, an `STR` tag and the repeat count: the curated count when a validation recorded one, otherwise `~` and ExpansionHunter's longer allele. A UPD reads `UPD chr2 · paternal`. The raw key, motif and per-member genotypes are in the tooltip. The cohort and standalone family pages now share one panel.

### Fixed
- Diagnostic All always showed a "Type" filter note, even with every type selected.

## [0.22.0] - 2026-10-20

### Changed
- **STR: flags follow the pipeline's new severity tiers** (ghfc-ngs `expansionhunter`). A genotype is flagged when its severity is high (pathogenic) or medium (premutation, or a pathogenic call of low confidence). Intermediate alleles, carriers and common risk alleles are info and no longer counted as flagged. This applies to:
  - the "EH (n flagged)" count;
  - the cohort flagged table and its download;
  - "common in cohort";
  - the REViewer message.

  Tables written before the pipeline had severities keep the old rule until their family is rebuilt.
- **Sample QC: contamination is the pipeline's own CHARR** (Lu et al. 2023). It is drawn against 0.02 (medium) and 0.05 (high). somalier's estimate is shown and exported as `somalier_charr`.

### Added
- **STR:**
  - `premutation` as a status and threshold (FMR1 55–200);
  - a Severity column on the family table and the cohort tables;
  - a `?` on low-confidence alleles;
  - Confidence and Read length columns in the locus dialog;
  - High and Med. counts on the cohort Loci table.
- **Sample QC: the new flags are worded.** Readers get neutral labels:
  - `CONSANGUINITY` and `RELATED_BEYOND_PEDIGREE`: "Pedigree relationship inconsistent";
  - `MZ_TWIN_CANDIDATE`: "Possible identical twins";
  - `UPD_CANDIDATE`, in a new "Uniparental disomy" category;
  - `RELATEDNESS_UNRELIABLE`, under "Across families".
- **Sample QC: pairs whose flags are only information** (twins, consanguinity, a third degree across families) are drawn as:
  - grey diamonds in the relatedness plots;
  - "Mismatch (info)" in the pair table;
  - dashed edges in the cross-family network, which gains a Flag column.

### Note
- Acknowledgements are keyed by flag code. Consanguinity acknowledged under `UNEXPECTED_RELATED_IN_FAMILY` comes back as an unacknowledged `CONSANGUINITY` (info) once the pipeline reruns.

## [0.21.0] - 2026-10-20

### Added
- **TSV downloads on every tab.** Every cohort tab except Statistics, and every family tab, has a Download button.
  - The file holds what the tab shows: the individuals that pass the cohort table's filter, or the family's selected members, plus the tab's own selectors. Every column is included, not just the one plotted.
  - Per-sample files start with `FID, IID, sex, phenotype`. A family file and a cohort file of the same kind have the same columns.
- **What each tab exports:**
  - **Ancestry:** region labels, PC1–PC20 and the admixture proportions (`snp_coverage`, `admix_<component>`, read from the pipeline's `Q.tsv`).
  - **Polygenic Scores:** the region and every trait z-score, one row per sample.
  - **ROH:** F_ROH, the percentile against the reference group, F short, intermediate and long, pattern, grade and flags. The family tab also offers the segments, the genes in ROH and kinship, separately or as a zip.
  - **STR:** the cohort's flagged genotypes; on the family page, every genotype of the selected members with validation and diagnostic.
  - **Sample QC:** three tables, separately or as a zip:
    - flags, with the curators' acknowledgement (status, note, user, date);
    - sex karyotype and sample quality;
    - pairs.
  - **Wombat and SVs:** the table's rows and visible columns as currently shown.
- **Files follow the user's role, like the tabs.** Readers get the neutral wording, an "Atypical" karyotype and only the pedigree-declared pairs. They get no evidence, observed relationship, ROH category or parents' kinship, and no acknowledgement note that could name a relationship.

### Fixed
- The Locus and Disease filters of the family STR table now filter.
- Filtering the family Polygenic Scores table keeps its sort.

## [0.20.0] - 2026-10-20

### Added
- **Sample QC from the ghfc-ngs `sampleqc` step (somalier).** The app reads each family's and cohort's `qc/<id>.<label>.{samples,pairs,flags}.tsv` and shows:
  - relatives across families;
  - relationship errors inside a family (not the declared parent, half vs full sibs, unexpectedly related, duplicates);
  - sex mismatches and possible sex chromosome aneuploidies;
  - contamination and Mendelian errors.
- **Family page: a Sample QC tab.**
  - The issues, with relatives in other families linked to their family.
  - A kinship matrix (observed against pedigree) and the KING kinship / IBS0-IBS2 plot.
  - X/Y copy numbers with karyotype boxes, and chrX heterozygosity.
  - Contamination and Mendelian error plots against the cohort.
  - A QC badge in the family header opens the tab, and flagged members get a warning icon.
- **Cohort page.**
  - A Sample QC column with a filter in the individuals table.
  - A Sample QC panel that follows the filters, with a link to somalier's HTML report.
  - A cross-family network where "hub" samples, related to many families, stand out. These are more likely contamination or population structure than real relatives.
- **Home cards** show how many individuals have a high-severity QC issue.
- **De novo warning.** Wombat de novo tables warn when the family has a non-parent, swap, duplicate, contamination or Mendelian flag.
- **Curator acknowledgements.** A curator can mark a flag as explained (known twins, known consanguinity…) or as a confirmed issue, with a note, in `validations/sampleqc.tsv`. Acknowledged flags stop counting in the badges, the cohort column and counts, the home cards and the de novo warning; confirmed issues still trigger the de novo warning.
- Wording follows the ROH tab: readers get neutral labels, no evidence and only the pairs the pedigree declares related.

## [0.19.2] - 2026-10-20

### Added
- **Suggested breakpoints in Read Evidence.** Each side of the SV gets:
  - **An exact position** from the carrier's soft-clipped and split reads, with a confidence: strong, moderate or weak. Reads whose orientation doesn't fit the SV type count less.
  - **The interval** the discordant pairs allow, shaded green.
  - **For deletions and duplications, a check against the coverage**: it must step the right way at that position. When the reads give no position (repeat-mediated events), the coverage step is offered as an approximate one, dashed and marked "≈". Inversions and other balanced events are judged on reads alone.
- **"Use as curated position"** fills the curated start and end from the suggestion ("Use approximate position" for a coverage step).
- **The read views follow the curated position.** Any change to it re-centres IGV, Read Evidence and the Copy Number View.
- **A 25 bp read-depth line** over the copy number in Read Evidence (2 × depth / the sample's median depth, not corrected by the panel) shows the exact point where coverage drops or rises.
- **"Show MAPQ 0 reads"** in Read Evidence draws evidence from reads with mapping quality 0 (repeats, such as the 22q11.2 low-copy repeats), faded. They never count toward the suggestions.

### Changed
- **WisecondorX calls get ±5 kb windows at each breakpoint** (was ±1.5 kb), since their breakpoints can be off by more than 1 kb.
- Split-read links rise higher, as half-ellipses with steep ends.

## [0.19.1] - 2026-10-20

### Changed
- **Read Evidence draws only the evidence**: discordant pairs (with both read bodies), split reads, soft-clip clusters and mates on another chromosome. Concordant pairs are no longer drawn; the range their insert sizes normally fall in (2.5-97.5 percentile) is shaded instead.
- **Copy number on the Read Evidence plot.** Each sample and window is one plot. Insert size is on the left axis, and copy number on the right, from the normalised coverage divided by the control panel, as in the Copy Number View. A sample with coverage but no CRAM keeps its copy-number row.
- **Split reads** join their two aligned parts with an arc rising above them, so the link no longer hides the deletion and duplication lines.
- **Soft-clip clusters** (several reads clipped at the same position, as at a breakpoint) are drawn as bars at the foot of the plot.
- **Read Evidence comes first in the SV dialog**, with the Auto / Single / Split toggle. The Copy Number View is folded and drawn when opened.
- **Navigation rulers** under each Read Evidence window. The mouse wheel no longer zooms Read Evidence or the Copy Number View: it scrolls the dialog, and dragging pans.

### Fixed
- The cohort that built a CNV control panel (GHFC-NT) now uses that panel itself, instead of showing "no control panel".
- The Read Evidence legend no longer wraps into the first plot, and the insert-size axis ends on a round value.

## [0.19.0] - 2026-10-19

### Added
- **Copy number in the SV dialog.** When the carrier has the normalised coverage bigWig written by the ghfc-ngs `cnv_coverage` step, the SV dialog shows a "Copy Number View" in place of the IGV coverage view. It has one plot per sample (carrier, parents, any added member) with a shared zoom.
  - Copy number = expected ploidy × sample / median of the cohort's control panel (`cnv_panel_dir` in its params file, built by the `cnv_panel` step), PAR-aware on chrX/chrY.
  - Bins where the panel is unreliable are greyed.
  - Each sample's karyotype and fit to the panel are shown, highlighted when the fit is worse than the panel's own cut-off.
  - Without a panel, or for a karyotype the panel has no group for, the sample's own normalisation is shown, with the reason.
  - Samples without a bigWig keep the IGV coverage view.
- **Read evidence (samplot-like) in the SV dialog.** Read pairs are drawn at their insert size and coloured by orientation, as in IGV's "pair orientation & insert size" colouring:
  - deletion (`+-`, large insert);
  - duplication (`-+`);
  - inversions (`++`, `--`);
  - mates on another chromosome.

  Split reads and a coverage strip are drawn too, for each sample with a CRAM. This needs the new `reference_fasta` config key, the FASTA the CRAMs were aligned to.
- **Single window or one per breakpoint.** The IGV read view and the read evidence show a short SV in one window, and a long one (wider than IGV's 30 kb alignment window) at each breakpoint. An Auto / Single / Split toggle overrides the choice.
- The sample dialog's coverage track uses the normalised bigWig when there is one.
- Workflow tab labels for the `melt`, `cnv_coverage`, `cnv_panel` and `sampleqc` steps.

### Fixed
- The read view of an SV shorter than 3 kb no longer shows two overlapping windows.

## [0.18.2] - 2026-10-19

### Fixed
- **Cohort de novo tables without a header** (CANDY_mpx, CANDY_mpx_ghent, PMS: a blank line, then the rows) are read with the header of a family file of the same config. This fixes:
  - the family Stats dialog's "Compared with the cohort", which showed every sample at 0;
  - the Search page, which built no table on them;
  - the cohort page's Statistics panel.

  With no family header that fits, each says the file has no header instead of showing wrong results.
- **The Search page's default Wombat source** no longer depends on file listing order, and a table under the older `.annotated.` name is listed under its config name (`strict`, not `annotated.strict`), in Search and in the Statistics panel.

## [0.18.1] - 2026-10-19

### Added
- **Variants per sample in the family Stats dialogs** (Wombat sub-tabs and SVs → WisecondorX), to spot a sample with an abnormal number of calls, such as a child with dozens of de novo variants.
  - "Variants per Sample (current view)": one bar per selected member, following the table's filters (SNV/Indel for Wombat, call type for SVs).
  - "Compared with the cohort (unfiltered)": each member's raw count placed on the cohort's per-sample distribution for the same config, from the cohort-level Wombat results or WisecondorX BED. Shows the group's median [Q1-Q3], the member's percentile and a High/Low flag (outside Q1-Q3 by more than 1.5 × IQR), with a histogram per group.
  - Founders and non-founders are compared within their own group, every pedigree sample counting 0 when it has no call, so the parents' zeros do not make every child with a de novo call look like an outlier.

## [0.18.0] - 2026-10-15

### Changed
- **Gene badges open GeneTrek's table on the gene** (`?search=SYMBOL`) instead of its gene dialog (`?gene=`).

### Added
- **More GeneTrek links.** The family ROH tab's "Genes in ROH" Gene column and the STR `Locus` columns (family EH table, cohort loci and flagged genotypes) are gene badges, coloured by your tiers and linked to GeneTrek; an STR locus is shown under its id and linked as its gene (`FXN_A` as FXN).
- **An SV's `chr:start-end`** (SVs tab, Search results) opens GeneTrek's table on the genes the region overlaps, through the `#search=` fragment so a patient's locus stays out of server logs.

## [0.17.0] - 2026-10-15

### Added
- **Gene badges come from GeneTrek.** The genesets and their membership are read from GeneTrek's API (`genetrek_url` in the config, default the dev instance `https://genetrek-py.dev.pasteur.cloud/`), no longer from the packaged `genesets.yaml` / `genesets.tsv`, which are removed.
  - The whole GeneTrek table is downloaded once per GeneTrek release (about 4 MB) and cached in `genetrek_cache_dir` (default `.genetrek/` next to the config). The app checks for a new release every 6 hours and works offline on the last one.
  - Old and alias symbols resolve to the current gene: a VEP `KIAA1109` is coloured as `BLTP1`.
  - Clicking a badge opens the gene in GeneTrek (`?gene=SYMBOL`). The tooltip lists the gene's genesets by tier, how many other GeneTrek sets it is in, and its LOEUF, s-het, EAGLE class and mode of inheritance.
  - The "+N genes" badge opens GeneTrek's table on all of the cell's genes, and the SV dialog links its region to GeneTrek. Both use the `#search=` fragment, which never reaches a server log.
- **Gene badge settings in the Profile page.** Each geneset is assigned to Tier 1-4, Show (tooltip only) or Hidden; a gene takes the colour of its best tier, and each tier's colour can be changed. A live preview shows the result. Settings are saved server-side in `<config name>.preferences.yaml`; an administrator can make theirs the default for everyone.
- **Geneset filters offer GeneTrek.** The Search page and the family Wombat tab filter on the user's tiers ("Tier 1", "Tiers 1-2"...), every GeneTrek gene list, and the data directory's `params/genesets/*.tsv` as before, in one grouped, searchable menu.
- **MOI column** in the Wombat tab, the Search results and the variant dialog: whether each variant fits its gene's modes of inheritance (✓ consistent, ? possible, ✗ inconsistent), with the reasoning in the tooltip.
  - Modes are pooled from all of GeneTrek's inheritance sources.
  - Recessive genes: homozygous with carrier parents, compound heterozygous in trans or in cis, a single heterozygous variant. Second hits are looked for among the rows of the same table.
  - A male's chrX call outside the pseudoautosomal regions is read as hemizygous, and a dominant variant inherited from an unaffected parent is marked possible.

### Changed
- Gene badges of every table are sorted best tier first and capped at six (the Wombat tab and validation pages were not), with a readable text colour on any badge colour.
- Gene badge tooltips in tables are multi-line instead of one clipped line.

### Fixed
- **Wombat sub-tabs opened more than 30 s after the page loaded showed no table.** Tables now wait for their container with no timeout, and each sub-tab is built the first time it is opened.
- The family Wombat tab's geneset filter split `VEP_SYMBOL` on `,` only, so `&`-joined symbols never matched.

## [0.16.0] - 2026-10-13

### Added
- **STR diagnostics**, beside the SNV and SV ones.
  - The family STR tab's locus dialog lists each member's diagnostic at the locus. A button opens the SNV/SV review dialog, with its history, ignore toggles, and a save form for curators.
  - Diagnostics are saved to `diagnostics/strs.tsv`, keyed `<chrom>:<pos>-<end>:<LOCUS>` like STR validations. `Gene` is the locus's gene and `Impact` the member's pipeline status.
  - A diagnosed member's cell in the EH table carries a `Dx` badge.
  - The Diagnostic All and Statistics pages list STR beside SNV and SV, with a type filter. So do the cohort and home summaries and the family Diagnostics panel.
  - Diagnostic kinds are listed once, in `utils/diagnostic_files.py`.

### Fixed
- **REViewer images were never found.** The ghfc-ngs `EXPANSIONHUNTER_REVIEWER` step writes them to the sample's `str/expansionhunter/reviewer/`, not `str/expansionhunter/`. It also names each image after the catalog's LocusId, so a variant of a multi-variant locus (`FXN_A`, `NOP56_CGCCTG`) is drawn under its locus (`FXN`, `NOP56`). The dialog now falls back to that name.
- **A missing REViewer image now says why.** The pipeline draws flagged genotypes only, so the dialog reports "not flagged, not drawn", "REViewer not run yet", or "no REViewer image" (an anomaly).
- **Diagnostic writes refuse viewers server-side**, not only by hiding the save button. The ignore toggle is disabled for them.
- **An SV diagnostic saved from the search page was never shown again.** The search page saved it under `chr:start-end`, while every table, the search page's own included, reads SV diagnostics under `chr:start-end:dup|del`. It is now saved with the type suffix, as the SVs tab does. Diagnostics already saved without the suffix are not migrated and stay hidden until rewritten.

## [0.15.0] - 2026-10-12

### Added
- **STR tab on the family page**, with an **EH** (ExpansionHunter) sub-tab, for the ghfc-ngs `expansionhunter` step. Other STR callers will add sibling sub-tabs.
  - **One row per locus, one column per member.** Each allele is a badge: its repeat count, coloured by its class against the STRipy thresholds, with `·F` for flanking-only support and `·IR` for in-repeat reads. Hovering a badge shows the CI, read counts, coverage and cohort percentile. Loci are sorted most severe first, and a switch hides loci where every member is normal.
  - **Flags:** `↑` marks a child whose longest allele is above both parents' at a non-normal status. **"common in cohort"** marks a locus where at least 5% of the cohort is pathogenic, carrier or intermediate (ARX_2, XYLT1), which points at the catalog or the thresholds rather than the person.
  - **Cohort context:** the cohort's `loci.tsv` and `genotypes.tsv` add the founders' median, 99.9th percentile and maximum, and each member's percentile. They are read off-thread and cached on the file's stat. Without them the tab shows the family alone.
  - **Locus dialog:** a histogram of the founders' longest allele with the thresholds and the members marked, then each member's full call. Each member's REViewer image is shown when the pipeline has written `<sample>.<catalog>.<locus>.svg`. The name pattern is set in `config/str.yaml`.
  - A family with several catalog labels gets a selector.
- **STR panel on the cohort page**, beside Statistics, Ancestry, Polygenic Scores and ROH, following the individuals table's filters.
  - **Loci:** how many of the individuals shown are pathogenic, carrier or intermediate at each locus. Beside them, the founders' 99.9th percentile and maximum and the "common" flag, both computed over the whole cohort.
  - **Flagged genotypes:** each shown individual's flagged calls, with allele badges and cohort percentile, linked to the family.
  - A locus opens the founders' distribution with the shown flagged individuals marked, the 25 longest at most.
- **STR validation for curators.** The family STR tab's locus dialog lists the locus's validations and, for curators and administrators, has a form: member, present / absent / uncertain, inheritance, curated repeat count and comment.
  - Validations are saved to `validations/strs.tsv`, keyed `<chrom>:<pos>-<end>:<LOCUS>`. The same file layout as `snvs.tsv` and `svs.tsv` means the shared loader and badges read it unchanged.
  - The writer refuses viewers server-side, not only by hiding the button.
  - A validated member's cell in the EH table carries the ✓ / ✗ / ? / ⚡ badge.
  - The All validations and Statistics pages list STR beside SNV and SV, with a type filter and their own colour. An STR row opens its family page.
- `ExpansionHunter (STR)` label for the step in the Workflow status. An STR flag in the sample dialog's availability badges.

### Fixed
- **The standalone family page (`/family/<FID>`) listed no members for a pipeline-written pedigree.** ghfc-ngs writes `<FID>.pedigree.tsv` without a header, and the loader took its first member for one, so the page reported "No pedigree file found for this family." and every tab had no members. Family pedigrees now use the cohort pedigree's header detection: headerless files are read by position (FID, IID, father, mother, sex, phenotype), headed ones by column name as before.
- **A headerless cohort pedigree whose first family ID starts with `FID`** (`FID001`) lost its first member to header detection, which was a prefix test. A first line is now a header only when it names a known pedigree column.

## [0.14.0] - 2026-10-12

### Added
- **A test suite** of 186 tests, running in about 5 seconds, plus 3 known issues marked `xfail`. It runs on a synthetic data directory and reference bundle built per test, so it needs no access to the data mount.
  - **Unit tests:** pedigree loading, cohort discovery and change detection; ancestry and ROH file discovery and readings; the run record and parameters file; curation file read/write; config loading; the plot builders; and validation of every shipped YAML config.
  - **UI tests:** NiceGUI's user simulation logs in as each role and opens the home, cohort and family pages and every lazy tab. A test fails on any error the page logs. The ROH tests check that relationship readings reach curators only.
  - **Boot test:** a subprocess starts the real server and checks `/health`, the login page and the redirect.
- **Continuous integration** (`.github/workflows/ci.yml`), on every pull request and push to main:
  - ruff check and format;
  - the tests on Python 3.12 and 3.13 with the locked dependencies;
  - a `package` job that builds the wheel, checks it ships its config and static files, installs it **without the lock file** (as `uvx` does in production) and runs the whole suite against it.
- **A weekly upstream check** (`upstream.yml`) does that unlocked install every Monday and opens an issue when an upstream release breaks the suite. The Polars 2.0 failure of 0.13.1 came from PyPI, not from a commit, and would have shown up here first.
- **Gated publishing:** `publish.yml` runs the full CI before building, and refuses a tag that does not match the `pyproject.toml` version.
- **Dependabot** for the workflow actions and `uv.lock`.
- `--show` option to open a browser tab on start, for local use.

### Fixed
- **The server no longer tries to open a browser when it starts.** `ui.run` defaulted to `show=True`, including in production on a headless VM.
- **Multi-line diagnostic comments corrupted `diagnostics/snvs.tsv`.** The comment was written raw, so a newline split the record. The SNV and SV validation dialogs had the same problem with tabs. All curation writers now share one `clean_tsv_field`.
- **A float-stringified parent ID (`0.0`, `-9.0`) was kept as a real parent**, so the founder was not counted as a founder. Sex and phenotype were already normalised.
- **`load_config` crashed on an empty key** such as `data_directories:`, and an empty `ancestry_reference_dir:` became the literal path `None`.
- **A hot-added data directory could reuse a removed directory's static prefix**, while the removed directory's route stayed mounted.
- **`norm_chrom("Chr1")` returned `CHR1`**: only a lowercase `chr` prefix was stripped.
- **A bundle release candidate (`1.0.0rc1`) outranked its release** when choosing the default variant.
- `pyyaml` is now a declared dependency. It was imported throughout but only arrived through NiceGUI.

### Changed
- **Polars 2 is supported again.** The requirement is now `polars>=1.17,<3`, and the lock moves to 2.0.0, which is what a production install now gets. 0.13.1 had capped it below 2 as a stopgap.
  - **The audit:** every Polars call was checked against the 2.0 upgrade guide, and the behaviours the app relies on were probed on 1.37 and 2.0 with identical results.
  - **Test runs:** the suite passes on 2.0.0, on the newest 1.x and on 1.17.
  - **Real data:** a pass under 2.0 covered every family and cohort tab and a cohort search.
- **Row order is now explicit wherever it matters.** Polars 2 guarantees no output order for a join or `unique` unless asked; the app had kept its order by accident until now. The two ancestry joins request `maintain_order="left"`, and the WisecondorX de-duplication in search requests `maintain_order=True`. `join(maintain_order=)` arrived in Polars 1.17, which is the new floor.
- **CI tests Polars 1.x alongside the locked 2.x:** a `polars-1x` job runs the newest 1.x and the 1.17 floor.
- **Plotly 7** (7.1.0) is now locked; the requirement is unchanged at `>=5`. **NiceGUI 3.18** and other minor updates come from Dependabot.
- **Ruff's rules are pinned** in `pyproject.toml` (`E4`, `E7`, `E9`, `F`, the set enforced until now). Ruff 0.16 widened its defaults, which would have changed what CI enforces on an upgrade. Markdown files are excluded from formatting. Ruff is locked at 0.17.
- **The Wombat aggregation** (one row per variant and sample) was duplicated between the family tab and search. It is now `utils/wombat.py::aggregate_by_variant`, and it joins distinct values in first-seen order. Before, the order of a merged cell such as `GENEA,GENEB` could differ between runs and Polars versions.
- The 13 pre-existing ruff errors are fixed, so lint can gate CI.

## [0.13.1] - 2026-10-12

### Fixed
- **Every cohort showed "Unreadable pedigree file — "column_6" not found" once polars 2.0 was installed.** Production installs with `uvx`, which ignores `uv.lock` and took polars 2.0.0 on its release. Polars 2.0 numbers headerless columns from `column_0` instead of `column_1`, and the pedigree loader renamed them by those generated names. It now renames by position, so it works on both versions.
- **The Wombat tab failed on every family under polars 2.0** ("`concat` was removed in version 2.0"). `str.concat` is replaced by `str.join`, available since polars 1.0.
- **The cohort Ancestry and Polygenic Scores tabs were silently empty for GHFC-NT and dataset_paper_20251010**, on any polars version. In those cohorts' tables, the first rows of `family_id` are numeric, so polars typed the column as an integer, then failed on `C0733-008-000-039` / `AU-ANG-FER-MPX3`. `read_tsv_or_none` takes a `string_columns` option, and the ancestry, PGS and ROH readers use it for their identifier columns.

### Changed
- `polars` is capped below 2 until the code is audited against 2.x. The fixes above make the known breakages version-independent, but a major release can also change behaviour silently.

## [0.13.0] - 2026-10-12

### Added
- **ROH tabs on the family and cohort pages** for the new ghfc-ngs `roh` step (runs of homozygosity, `ancestry-pgs roh` ≥ 0.2.0). Every value is shown against the unrelated 1kGP+HGDP reference genomes of the person's own ancestry region, read from the bundle's `roh/` tier (`v<bundle>/bundle/roh/reference_froh.tsv.gz`), before the fixed published thresholds. Region comes first because the fixed thresholds alone would mislabel whole populations: among unrelated reference genomes, F_ROH > 0.03 holds for 150 of 682 South Asian and 61 of 126 Middle Eastern ones, against 17 of 665 European ones.
- **Family ROH tab**, with five blocks:
  - **Summary:** one row per member, with F_ROH and its percentile within the region, F by segment class, the population-vs-recent pattern, a grade, and flags for a UPD candidate chromosome or a failed call-rate check.
  - **Against the reference:** per reference group, a strip of the reference with the members on it (F_ROH, F long or F short + intermediate), and the population-vs-recent plane whose 95th-percentile lines delimit the pipeline's four patterns.
  - **Genome:** chromosomes 1–22 with each member's segments by length class and the reference hotspots shaded; a UPD candidate chromosome is outlined.
  - **Genes in ROH:** longest ROH first, with genes mostly inside reference hotspots hidden by default.
  - **Kinship:** KING kinship checked against the pedigree, so a declared parent–offspring or sib pair at the wrong degree is flagged as a mismatch.
- **Cohort ROH sub-tab** in the right-hand panel, following the individuals-table filters:
  - horizontal split violins per reference group, the cohort beside its reference, or the cohort alone by phenotype or sex;
  - the population-vs-recent plane for one group;
  - grade and pattern counts;
  - a list of flagged individuals linked to their family;
  - a full-width dialog with the founders' genome-wide ROH islands (`roh_windows.tsv`) against the reference hotspots.

  Only the cohort `froh` table is read: the cohort segment and gene tables run to 31 and 93 MB.
- **Role-dependent wording.** Viewers see neutral grades (Elevated / High / Very high / Mild) and only the kinship pairs the pedigree declares related. Curators and administrators also see the relationship each grade is typical of, the parents' kinship and expected F, and every kinship pair. A first-degree reading can reveal incest, which ACMG asks to be handled under an institutional policy. This is presentation only: the raw TSVs remain reachable through the data directory's static route.
- **`config/roh.yaml`** — grade thresholds and wording, pattern and segment-class colours, KING degree bounds. They mirror ancestry-pgs constants; no component hard-codes a threshold.
- `roh` in `config/pipeline_steps.yaml`, labelled "Runs of homozygosity".

### Changed
- Ancestry file discovery takes the set of result kinds as a parameter, so the `roh/` directories, which share the `<id>.apgs_b<bundle>_<tag>.<kind>.tsv` naming, reuse it.
- The chromosome backbone of the search statistics ideogram (bands, outline, Mb grid) moved into `utils/ideogram.py`, shared with the ROH genome view. Its output is unchanged.
- The phenotype and sex plot palettes moved from the cohort PGS tab into `utils/pedigree_labels.py`.

## [0.12.0] - 2026-09-22

### Added
- **A cohort is now any directory under `cohorts/`.** Directories whose `<NAME>.pedigree.tsv` is missing — or present but unparseable — are listed on the home page as *Incomplete Cohorts*: amber, non-clickable cards naming the file they need, instead of being silently skipped. Six such directories existed on the GHFC data directory and were invisible. `DataStore.cohorts` still holds only explorable cohorts, so the header dropdown, validation pages and search are unchanged; the rest arrive in `DataStore.incomplete_cohorts`, which distinguishes a missing pedigree from an unreadable one because the two need different fixes.
- **Full-width `General` / `Workflow` tabs on the cohort page.** `General` is the individuals table and the right-hand panels that follow its filters, unchanged. `Workflow` gathers everything about the pipeline run and exists only for a cohort with a `<NAME>.params.yml`.
- **Pipeline status, read from the workflow's own run record.** The ghfc-ngs workflow writes `.ghfc-ngs.state.json` into each cohort directory with per-step `{done, total, pct}` counts (schema: `COHORT_STATE.md` in bourgeron-lab/ghfc-ngs), and genetics-viz never recomputes them — re-deriving would duplicate the pipeline's planning logic, and a cold existence check costs ~5.6 ms over the CIFS mount, so a 2700-sample cohort would be a minute of stat calls. Home page cards gain a foldable **Status** item: a status chip and an "N/M steps complete" summary collapsed, per-step progress bars when opened. The Workflow tab shows the full record — run timings and duration, last successful run, pipeline version and commit, run history — with warnings for a stale `running` record (there is no heartbeat, so a killed run leaves one behind), counts measured before the run, and outputs a failed run may not have published. A step the record reports as `null` renders "not measured", never 0%: unmeasured and zero are different facts.
- **Drift detection.** The Workflow tab re-hashes the local pedigree and params file against the checksums in the record, flagging "the pedigree has been edited since this run". The paths inside the record are compute-cluster paths, so the local files are hashed instead.
- **The cohort's `params.yml` beside the pipeline's own reference for it**, two side-by-side panels under the status, each scrolling internally so the page never grows a long scrollbar. The parameters file is shown verbatim — these files carry a great deal of their meaning in comments. The reference is `documentation/params.md`, fetched at runtime from the public ghfc-ngs repository (stdlib `urllib`, 10s timeout, cached an hour — no new dependency) rather than vendored, since it changes often. Read outside its own repository, its repo-relative links are rewritten to absolute GitHub URLs opening in a new tab, with `/tree/main/` for the one directory link; its 73 in-page anchors resolve inside the panel via markdown2's `header-ids`, with clicks intercepted so the page neither jumps nor collects a URL fragment. A failed fetch degrades to a note naming the error plus a link to the original, leaving the rest of the tab intact.
- **`config/pipeline_steps.yaml`** — step labels, display order and run-status chip colours. Display metadata only: completion is measured by the pipeline, never recomputed here. A step the pipeline reports but this file does not know is shown under its raw name.

### Fixed
- **An unknown cohort name in the URL returned a 500.** `create_header` passed the name straight into a `ui.select` value; Quasar raises when the value is not among the options, which fired before the pages' own "cohort not found" branches could run, making them dead code. Those branches now work, and a pedigree-less cohort directory explains what it is missing instead of erroring.
- **The cohort cards lost most of their width to layout defects**, which the new Status panel made visible: `nicegui-card`'s own 16px padding doubled with the 16px on each `q-card__section` (32px a side), and because `nicegui-card` is a flex column with `align-items: flex-start` the sections did not stretch, coming out 274px wide in a 320px card — so nothing inside lined up with the card's edges. Sections now carry `w-full`, the card carries `p-0`, and padding is `pl-2 pr-3` to offset the 4px accent border. Usable content went from 250px to 296px with symmetric insets, and no element overflows at any container width from 360px up.
- **Two flex-sizing overflows on the new panels.** Long file paths escaped their card because `.nicegui-column` sets `align-items: flex-start`, so a child is sized by its content and `truncate` has no width to act on; those labels now carry `w-full`. Separately, the step rows' `whitespace-nowrap` counts gave the home card's Status expansion a 299px min-content width — a flex item will not shrink below its content, and that floor propagates up the ancestor chain — which forced the card section to 331px inside a fixed 320px card and clipped the phenotype table's last column. Every flex item in the chain now carries `min-w-0` and truncates, and the card's compact rendering drops the unit noun (`100% · 522/522` rather than `100% · 522/522 samples`); the Workflow tab, which has the room, keeps it.

### Changed
- Cohort discovery was hand-copied into `DataStore.load()`, `take_snapshot()` and `reload()`; all three now share one `_iter_cohort_dirs()` generator. `take_snapshot()` additionally watches `<NAME>.params.yml` and keys every cohort directory rather than only pedigree-bearing ones, so a directory appearing or gaining a pedigree is detected. The run record is deliberately *not* watched: the pipeline rewrites it on every run, and a reload re-parses every pedigree in the data directory — its reader caches on the file's own stat instead.

## [0.11.0] - 2026-09-21

### Added
- **Ancestry and Polygenic Scores tabs on the cohort page** — the right-hand Statistics panel becomes three tabs. Both new tabs read the cohort-level pipeline output and follow the filters applied to the individuals table.
- **Cohort Ancestry tab** — every cohort member on PC1×PC2 over the reference panel, which a checkbox (ticked by default) can hide. Points are coloured by predicted region and outlined red when flagged as outliers; the tooltip names the sample, region, population and confidence. Beside the plot, the region distribution as a pie plus a count/percent table. For EAGER that is 586 samples across 7 regions.
- **Cohort Polygenic Scores tab** — violins of the cohort's z-scores for a selected trait. The grouping axis is Phenotype, Sex or predicted Region; an optional split draws each group as two sex halves; and a region multiselect restricts which samples are included. Group sizes are printed on the axis and the individual points are overlaid, because filtering by region can leave groups of a handful of samples.
- **Bundle/filter-tag selector** on the cohort and family ancestry tabs, shown when a directory holds more than one variant.

### Fixed
- **The strictest filter-tag variant is now preferred.** Ancestry directories now hold two variants of the same bundle — `dp6gq15` and `dp10gq20` — and the tags were sorted as plain strings, so `dp6gq15` won because "6" sorts above "1" at the third character. The family tabs released in 0.10.0 therefore switched variant silently as soon as the second one appeared. Tags are now ranked by their depth and GQ thresholds.
- Pedigree sex and phenotype codes written float-stringified (`-9.0`, `1.0`, as the PMS cohort does) were not recognised as their integer equivalents, so they grouped as unknown.

### Changed
- The shared PCA plot builders moved from the family Ancestry tab into `utils/ancestry_plots.py`, and the phenotype label helpers from `home.py` into `utils/pedigree_labels.py`, so the cohort tabs reuse them rather than importing another page's private functions.
- The cohort violins are Plotly, the second deliberate exception to the project's ECharts default — ECharts has no violin series, and the sibling Statistics tab is already Plotly.

## [0.10.0] - 2026-09-17

### Added
- **Ancestry tab on the family page** — the family's samples projected onto the versioned `ancestry-pgs` reference panel. Four scatter plots: PC1×PC2 and PC3×PC4 over the full panel, then the same pair zoomed to a window covering the family plus every reference sample of the regions predicted for it (~13× magnification for a single-region family). The reference cloud is one series per region, so the legend toggles regions, and it is `silent` so it never steals a tooltip from a family point. Family samples are labelled diamonds filled with their predicted region colour, red-bordered when flagged as outliers. Above the plots, a table of predicted region, population, both confidences, kNN distances, local reference density and the outlier flag.
- **Polygenic Scores tab** — the per-family `pgs_zscore.tsv`, which the pipeline writes wide (one row per sample, ~100 score columns), transposed to one row per trait and one column per member, with a diverging colour scale on the z-scores and a text filter on the trait name.
- **`ancestry_reference_dir` config key** — root of the references tree. The bundle version subdirectory is derived per family from its own filename tag (`apgs_b1.0.0` → `v1.0.0`), so each family is plotted against the panel it was actually projected onto rather than a globally pinned one. The reference panel is cached per bundle and read once, not per page load.
- **`config/ancestry.yaml`** — colourblind-safe reference region colours and the PGS z-score colour scale.
- Both tabs carry existence probes and are disabled for families with no `ancestry/` directory, which is currently the majority. They are wired into the cohort-scoped and the standalone family page.

### Fixed
- **`save_config` no longer drops unknown top-level keys.** It rewrites the whole YAML from known fields, so `poll_interval` was already being silently deleted the next time an administrator added a user or a data directory; `ancestry_reference_dir` would have inherited the same fate. Both are now emitted when they carry information.
- A missing or unconfigured reference panel is treated as a degraded state rather than a failure: the Ancestry tab still plots the family's own samples and names the missing config key, instead of erroring.

### Changed
- CLAUDE.md structure listing corrected — the continuous-score config file is `config/continuous_scores.yaml`, not `config/score_colors.yaml`.

## [0.9.2] - 2026-09-16

### Added
- **Tabbed family information panel** — Notes and Diagnostics moved into a **General** tab, joined by **Ancestry** and **Polygenic Scores** placeholders. General opens by default, and a half-typed note survives a visit to another tab.
- **`/health` endpoint** — unauthenticated liveness probe reporting the running version from `importlib.metadata`, so a deployment can assert that the release it asked for is the one answering. Left unauthenticated on purpose so a pool monitor can use it instead of a bare TCP connect.

### Changed
- **Pedigree table is much shorter** — row pitch dropped from 73px to 29px and the table from 640px to 260px for an eight-member family, by removing the inter-row flex gap and using a dense checkbox.
- The member table became a shared `render_member_selector` component, so the standalone family page — which carried a near-verbatim copy with the same bugs — is fixed by the same change.

### Fixed
- **Pedigree table alignment** — the grid style was applied to the header row and each member row separately, making every row an independent grid sized to its own content. Six of seven columns drifted by up to 75px. A single grid container now holds every cell.
- **Diagnostics merged across curators** — a variant recorded by two curators appeared twice; it is now one row, with a tooltip listing each curator's own verdict and date.
- **Diagnostics merged across individuals** — a variant reaching the same verdict in several family members is one row listing every barcode in pedigree order. Members that reached *different* verdicts stay on separate rows.
- **Empty wombat files no longer log as failures** — wombat leaves a zero-byte file when a step had nothing to report (a family without two sequenced parents has no de novo mutations) and polars raises `NoDataError` on it, which was being caught by the handler meant for real failures and logging a traceback on every visit. A new `read_tsv_or_none` helper returns `None` for such a file, leaving `FileNotFoundError` to propagate.

## [0.9.1] - 2026-08-25

### Fixed
- **Pinned TanStack table-core to 8.21.3 on the CDN** — the jsDelivr URL was unpinned, so it resolved to whatever npm tagged `latest`. TanStack Table v9.1.2 became `latest` and dropped the `getCoreRowModel` export, so every table threw `TanStack.getCoreRowModel is not a function` at load time and never rendered.

## [0.9.0] - 2026-05-20

### Added
- **Top-level tab disabling** — the Wombat and SVs tabs on the family page are greyed out and unclickable when no data files exist for a family.
- **Sub-tab disabling** — within Wombat, sub-tabs whose DataFrame is empty are disabled rather than showing blank content, distinguishing "step ran but found nothing" from "step was not run".
- **Variant count badges** — enabled sub-tabs show a count, e.g. `config_name (42)`, so the volume is visible without opening the tab.
- **Smart default selection** — the first tab with actual data is auto-selected, skipping disabled ones.
- Lightweight `probe_wombat_data()` / `probe_svs_data()` helpers check for data existence without a full DataFrame load, so tab state can be decided before rendering.

## [0.8.1] - 2026-05-19

### Fixed
- **Deferred config save until the password dialog is closed** — in admin user management, `save_config()` ran before showing the generated password, so the app could reload before the administrator had copied it. The YAML write now happens only when the dialog is dismissed.

## [0.8.0] - 2026-05-19

### Added
- **Live data refresh** — filesystem changes are detected automatically; no restart when cohorts are added, removed or updated on disk.
- **Background polling** — one server-wide task checks all data directories every 30 seconds for new, removed or modified cohort pedigree files, using `asyncio.to_thread()` for the stat calls so the event loop stays responsive.
- **Toast notifications** — connected browsers are told what changed.
- **Manual refresh button** in the header for an immediate check.
- **Atomic reload** — `DataStore` rebuilds its cohort dict in a local variable before swapping, so concurrent page loads always see a consistent state.
- **`poll_interval: <seconds>`** config key to adjust the polling interval (default 30).

## [0.7.9] - 2026-04-08

### Added
- **Help menu in header bar** — small `help_outline` icon button (right side, before the user menu) with three options: Email (mailto to the Teams channel address with "Bug report" subject), Teams channel (opens the genetics-viz channel in a new tab), and Report bug (opens a new GitHub issue in a new tab).

## [0.7.8] - 2026-04-08

### Fixed
- **Family page member table** — rebuilt the member selection table using native NiceGUI CSS-grid layout instead of a static HTML table + JS DOM moves. The previous approach was fragile because Vue/NiceGUI reconciliation would revert the manual moves. Checkboxes and "only" buttons now render reliably inside their cells. Same fix applied to standalone family page.

## [0.7.7] - 2026-04-08

### Fixed
- **Family page member table checkboxes** — fixed timing bug where checkboxes and "only" buttons sometimes rendered below the table instead of inside the Select column. The DOM-move JavaScript now polls for both source elements and target cells with bounded retries.
- **DataTable container race condition** — replaced fragile MutationObserver with a bounded polling loop in `tanstack_table.py`. Fixes intermittent issue where the Wombat/SVs dataframe wouldn't render and required a page reload.

## [0.7.6] - 2026-04-07

### Changed
- **Cohort card phenotype table styling** — shortened phenotype labels (`2 (aff)`, `1 (unaff)`, `-9 (unk)`) with `whitespace-nowrap` to prevent wrapping; wrapped the table in a rounded border with a separator line between the header row and data rows.

## [0.7.5] - 2026-04-07

### Changed
- **Cohort card phenotype breakdown** — replaced the single Affected/Diagnosed columns with a per-phenotype table showing N samples, Pat/Unc counts, and diagnostic yield % for each phenotype category (2=affected, 1=unaffected, -9=unknown). Only categories with samples are displayed.

## [0.7.4] - 2026-04-07

### Changed
- **Cohort card diagnosed count** — now shows two numbers: pathogenic (red) / uncertain-only (amber), where uncertain-only counts samples with at least one uncertain diagnostic but no pathogenic. Percentage still based on pathogenic / affected.

## [0.7.3] - 2026-04-07

### Added
- **Family navigation arrows** — prev/next arrow buttons on the cohort family page header to navigate between families in the current project. Hidden for single-family cohorts; previous hidden on first family, next hidden on last.

## [0.7.2] - 2026-04-07

### Added
- **Family Notes panel** — new panel above Diagnostics on both cohort and standalone family pages for adding free-text notes to families or specific samples. Notes stored in `notes/notes.tsv` with hard delete support. Includes add form (with optional sample selector), delete button, info tooltip, and write protection.

## [0.7.1] - 2026-04-03

### Added
- **Version number in header** — displays version below the app name in the header bar
- **Cohort card statistics** — home page cohort cards now show affected count, pathogenic diagnostic count, and diagnostic yield percentage alongside families and samples
- **Diagnostic column in cohort table** — new Diagnostic column with colored badges (pathogenic/uncertain/benign) and multiselect filter (including NA for undiagnosed) on the cohort detail page

## [0.7.0] - 2026-04-03

### Added
- **Quick Search on home page** — search samples by barcode or families by FID with async O(1) sharded filesystem lookup; results show data availability badges and action buttons
- **Sample visualization dialog** — fullscreen IGV.js viewer with bedgraph (CNV), CRAM (alignments), and VAF tracks for any sample, with locus navigation input
- **Standalone family page** (`/family/{fid}`) — browse family data independently of any cohort/project, with data availability panel, member table, diagnostics, and analysis tabs (Wombat, SVs)
- **Standalone pedigree parser** (`utils/pedigree.py`) — parses per-family pedigree files with extended column name support (FatherBarcode, MotherBarcode, Pheno_*)
- **Data availability checker** (`utils/data_availability.py`) — checks existence of data files for samples and families

### Changed
- `show_sv_dialog()` and `show_variant_dialog()` accept optional `family_members_override` and `sample_parents_override` parameters for standalone pages without cohort context (backward compatible)
- **Sharding detection** — replaced "all children must be single-char" heuristic with try-sharded-first strategy to support hybrid directories with both shard buckets and direct entity folders

## [0.6.9] - 2026-03-24

### Added
- **CRAM split view re-centering** — when ROI is updated (via suggestion click, curated position refresh, or New Start/End buttons), the Read-Level Split View automatically re-centers both panes on the new boundary positions (±1500bp windows)

## [0.6.8] - 2026-03-24

### Fixed
- **IGV.js ROI display** — reverted from v3.8.0 to v2.15.13 (latest 2.x) due to ROI rendering incompatibility in v3. ROI overlays now display correctly again.

## [0.6.7] - 2026-03-24

### Fixed
- **IGV.js v3 ROI display** — ROI regions were not rendering because `loadROI()` became async in v3. Moved ROI loading from `createBrowser` config to explicit `await browser.loadROI()` calls after browser creation. Also wrapped dynamic `_update_roi()` in an async IIFE.

## [0.6.6] - 2026-03-24

### Added
- **SV suggestion inheritance preselection** — clicking a coordinate suggestion from a parent now auto-selects "paternal"/"maternal" in the inheritance dropdown
- **Dual overlap percentages** — SV suggestions now show both "ours" (overlap as % of current SV) and "theirs" (overlap as % of suggested SV) for at-a-glance size comparison

### Changed
- **IGV.js upgraded from v2.15.11 to v3.8.0** across all 6 page files (7 script tags)

## [0.6.5] - 2026-03-23

### Fixed
- **SV ideogram colors in stats dialog** — SVs from to_validate files with non-standard `call` column values (e.g. "primary") now correctly show GAIN/LOSS colors on the ideogram. Uses shared `infer_sv_type()` for consistent SV type detection across all code paths.

## [0.6.4] - 2026-03-23

### Added
- **SV coordinate suggestions** — panel between SV Details and CNV Coverage View suggesting curated coordinates from overlapping validated SVs; prioritizes parents, then family, then cohort; click to auto-fill curated position fields and update ROI

## [0.6.3] - 2026-03-23

### Fixed
- **SV validation display for GAIN variants** — validations on "gain" type SVs now display correctly in the validation file table. Root cause: inconsistent SV type inference between save (sv_dialog.py, checked only `call` column) and lookup (file.py, checked `wisecondorX` first). Extracted shared `infer_sv_type()` utility in `utils/wisecondorx.py` used by all callers.

## [0.6.2] - 2026-03-16

### Added
- **VAF View panel** in SV dialog — displays VAF bedgraph tracks (`*.vaf.bedgraph.gz`) as a collapsible panel between CNV Coverage View and Read-Level Split View, rendered as scatter points with fixed 0–1 Y-axis scale

## [0.6.1] - 2026-03-16

### Added
- **Sharded directory support** — two-level sharding for `samples/` and `families/` directories (`utils/sharding.py`), with auto-detection and backward compatibility for flat layouts

### Changed
- All sample/family path constructions now use `get_sample_path()`/`get_family_path()` and `get_sample_url()`/`get_family_url()` from the sharding module

## [0.6.0] - 2026-03-09

### Added
- **Authentication & authorization** — YAML-configured user accounts with SHA-512 password hashing and role-based access control (reader, curator, administrator)
- **Login page** — username/password login with session persistence
- **Profile page** — view current role and change password
- **Admin: user management** — add, remove, change role, and reset password from the web interface
- **Admin: data directory management** — add, remove, and set default data directory from the web interface
- **Multi-data-directory support** — per-user data directory selection via header dropdown; YAML config lists multiple directories with descriptions and default flag
- **YAML config file** — single config file replaces CLI data_dir argument; holds data directories, user list, and auto-generated storage secret
- **Write protection** — save operations (validations, diagnostics) gated by `can_write()` at both UI and backend level
- **Config model** (`config_model.py`) — dataclasses for config loading/saving with file locking and password helpers

### Changed
- **CLI** — accepts `config_file` (path to YAML) instead of `data_dir`
- **App init** — loads config, initializes multi-store registry, registers per-directory static files, sets storage secret from config
- **Header** — data directory dropdown replaces reload button; user menu with profile/logout added
- **Static file URLs** — migrated from hardcoded `/data/` to `get_static_prefix()` for per-directory routing
- **All page handlers** — auth guard added (`check_auth()` redirect at top of every page)
- **Validation/diagnostic saves** — username read from session instead of OS user; `can_write()` guards added

## [0.5.0] - 2026-02-24

### Added
- **Exclude samples filter** — search Individuals panel: enter sample IDs (space/comma/blur) as deletable chips to exclude from results
- **SV support in statistics dialog** — SV deduplication, gain/loss classification in consequence chart, SV checkbox filter, variant type breakdown in subtitle
- **SV validation file support** — auto-detects SV format in validation files, loads `svs.tsv`, opens SV dialog instead of IGV
- **SV type inference** — extracts SV type from wisecondorX call, generic call, type column, or ratio sign
- **WisecondorX tooltips** — CNV call cells show ratio/zscore on hover in search results
- **Numeric column formatting** — integer columns display with thin-space (U+202F) thousands separator and proper Unicode minus sign
- **Column type config** — `type: int/float` in `column_names.yaml` auto-infers numerical sorting and number cell rendering

### Changed
- Sex and Phenotype dropdowns now on same row in search Individuals panel
- Validation file table preserves sorting and pagination across refreshes
- Statistics dialog subtitle shows variant type breakdown (SNVs, Indels, SVs)

### Fixed
- `KeyError: 'error'` on ratio Min/Max input fields — `props(remove="error")` was deleting the key from internal props dict

## [0.4.1] - 2026-02-18

### Added
- **Inheritance inference** — variant dialog auto-infers inheritance mode from VCF genotypes
- **Search TSV export** — download search results as TSV or save as validation file
- **Validation file stats** — stats dialog (chromosome distribution, consequence/status charts, ideogram) now available on validation file pages
- **Validation file badges** — gene, consequence, ClinVar, and score badges in validation file tables
- **Validation file column selector** — column visibility dialog with presets for validation pages
- **CRAM version check script** — `utils/check_cram_version.sh` checks IGV.js compatibility via magic bytes
- **FID column** added to all view presets

### Refactored
- Extracted `utils/locus.py` — centralized locus query parsing and filtering (from search.py)
- Extracted `utils/genesets.py` — shared geneset loading (was duplicated in 3 places)
- Extracted `utils/wisecondorx.py` — WisecondorX parsing, CNV classification, and color utilities (from svs_tab.py)
- Extracted `components/search_stats.py` — reusable stats dialog (from search.py)

### Changed
- Validation all/statistics pages now async with offloaded I/O

## [0.3.0] - 2026-02-16

### Added
- **TanStack Table** - New high-performance DataTable component replacing NiceGUI ag-grid, with virtual scrolling, column sorting, and custom cell renderers
- **Cohort-wide Search** - Tabbed search parameters panel with Variants and Individuals tabs
  - Individual filters: sex, phenotype, and "only samples with both parents"
  - Filters use pedigree data from the Cohort object (no duplicate parsing)
- **Variant Statistics dialog** - Stats button on search results showing:
  - Stacked bar chart of variants per chromosome by validation status
  - Consequence distribution pie chart (highest priority per variant)
  - Validation status distribution pie chart
  - Interactive SVG ideogram with cytoband rendering and variant positions
  - SNV/Indel type filter checkboxes with live refresh
- **ClinVar `Likely_risk_allele`** term added to clinvar_colors.yaml
- **Column names config** (`column_names.yaml`) for display labels, groups, sorting, width, and drop flags
- **View presets** for quick column visibility switching
- **Gene scoring** utilities and score color coding
- **Stats panel** component for carrier frequency box/bar plots
- **CLAUDE.md** project-level instructions for Claude Code

### Changed
- **Pedigree parsing** now handles `-9` as missing/unknown (in addition to `0` and empty)
- **Pedigree header detection** strips leading `#` (supports `#FID` headers)
- **Search panel** uses dense Quasar props for a more compact layout
- Pedigree data loaded from the already-parsed `Cohort` object via `DataStore` instead of re-reading the file

### Refactored
- Extracted `utils/vep.py` - consolidated VEP consequence utilities from 3 files
- Extracted `utils/clinvar.py` - consolidated ClinVar utilities from 2 files
- Extracted `utils/view_presets.py` - view preset logic from wombat_tab
- Extracted `utils/cytobands.py` - cytoband/ideogram constants shared between wombat_tab and search
- Fixed `header.py` coupling - no longer imports from page components

### Fixed
- **TanStack table memory leak** - replaced global `_injected_pages` set with per-client attribute
- **MutationObserver timeout** - 30-second timeout prevents indefinite observation
- `Sample.is_founder` now correctly handles `-9` parent IDs

### Removed
- Dead code: `variant.py` (unused variant page)
- Dead code: `waves_backup.py` (obsolete backup)

## [0.2.0] - 2026-01-20

### Changed
- **Breaking Change**: New validation TSV format with `Comment` and `Ignore` columns
  - Use migration script: `utils/snvs_validations_migration_0.1.1_to_0.2.0.sh`
- Validation form now defaults to "present" status
- Validation/all page aggregates by Variant/Sample with unique user lists
- "in phase MNV" normalized to "present" for conflict detection but displayed distinctly

### Added
- Enhanced inheritance options: "not maternal", "not paternal", "homozygous"
- New validation status: "in phase MNV" for multi-nucleotide variants
- Comment field for free-text notes on validations
- Ignore functionality to exclude specific validations from statistics and conflict detection
- Interactive validation guide accessible via info button in dialog
- Validation history shows ignore toggle for each entry
- Ignored validations displayed with reduced opacity
- Statistics page excludes ignored validations and shows separate ignored count

### Fixed
- Fixed context issues when toggling ignore status from within dialog
- Fixed table refresh after validation changes

## [0.1.1] - 2026-01-15

### Added
- Initial stable release
- Core validation functionality
- IGV.js integration
- WAVES validation support
- Multi-cohort management
- Family structure visualization
- DNM and WOMBAT analysis tables
- Variant validation tracking with inheritance patterns
