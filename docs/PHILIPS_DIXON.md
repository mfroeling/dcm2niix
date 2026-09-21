# Philips Dixon Water/Fat/In-Phase/Out-of-Phase Handling

Read this before touching Philips Enhanced MR volume splitting, dimension-index
handling, or anything that looks like a "which type is this frame" decision.

## Coverage tested

Fixes here were verified against real Philips Enhanced MR data spanning
software versions 5.1.x/5.7.x through 12.3.0, covering: classic and Enhanced
MR SOP classes; Dixon source series (Magnitude/Real/Imaginary/Phase echoes)
and Dixon recon series (Water/Fat/In-Phase/Out-of-Phase); plain multi-echo
gradient-echo and multi-echo spin-echo (T2 mapping) series with no Dixon
option set; and DTI/DWI series. Also verified clean against the project's
in-tree regression suite (`dcm_qa`, `dcm_qa_nih`, `dcm_qa_uih`).

## The problem

Philips Dixon acquisitions (mDIXON / GRE-Dixon) produce two kinds of series:

- **Source**: the raw acquired echoes, as Magnitude/Real/Imaginary/Phase —
  dcm2niix already understood this vocabulary before this change.
- **Recon**: derived Water, Fat, In-Phase and Out-of-Phase maps computed
  on-scanner from the source echoes.

dcm2niix had no concept of the second group. `ComplexImageComponent` reports
`MAGNITUDE` for all four recon types, so they were indistinguishable from each
other and from an ordinary magnitude image using the tags already parsed.
Symptoms:

- Water and Fat (same TE, no other distinguishing tag) silently merged into
  one multi-volume file instead of two separate outputs.
- In-Phase and Out-of-Phase either merged, or one clobbered the other's slot,
  depending on which frames happened to share an echo/type bucket.
- Output filenames and JSON sidecars carried zero information about which
  volume was which — no `_real`/`_imaginary`/`_ph`-style marker existed for
  these four types, so downstream tools could only guess from echo time.

## Where the real type lives

Enhanced MR Philips files carry the per-frame type in
`MRImageFrameTypeSequence` (0018,9226) → `FrameType` (0008,9007), last array
element. The vocabulary changed between software releases:

- Pre-R12: short codes `W`, `F`, `IP`, `OP`.
- R12+: long codes `WATER`, `FAT`, `IN_PHASE`, `OUT_OF_PHASE`.

Both eras *also* duplicate this as a short code in the legacy per-frame
compatibility block `(2005,140F)` → nested `(0008,0008) ImageType`, e.g.
`DERIVED\PRIMARY\W\W\DERIVED`, and in the private tag `(2005,1011)`. This
legacy block is present in every file checked across the full version range,
and is already walked by the existing classic `ImageType` parser (the same
underscore-joined string matching used for `_R_`/`_M_`/`_I_`/`_P_`). That is
the anchor this fix uses: it is the one representation that has not changed
across a decade of Philips software, so no per-version vocabulary table is
needed — only the short-code tokens `_W_`/`_F_`/`_IP_`/`_OP_`.

## Fix 1 — recognize the tokens (`// start/end dixon label fix`)

Mirrors the existing `isReal`/`isImaginary`/`isPhase`/`isMagnitude` machinery
one-for-one, gated behind `d.manufacturer == kMANUFACTURER_PHILIPS`:

- `nii_dicom.h`: four new `TDICOMdata` booleans (`isHasWater`, `isHasFat`,
  `isHasInPhase`, `isHasOutPhase`), four new `TDTI4D` per-frame arrays.
- `nii_dicom.cpp`: token detection in the `ImageType` substring block, reset
  alongside the others when `ComplexImageComponent` is (re-)parsed per frame,
  and — critically — fed into the per-frame `imageType` index used to keep
  volumes apart (values 4–7, after the existing 0=magnitude/1=real/2=imaginary/3=phase).
- `nii_dicom_batch.cpp`: `isSameSet()` (classic per-instance stacking) and the
  `gradDynVol` equality test inside `saveDcm2Nii()` (intra-file volume
  splitting for Enhanced MR) both gained a check on the four new flags, so
  two frames that differ only by Dixon type are no longer treated as repeats
  of the same volume.
- New human-readable filename suffixes `_water`/`_fat`/`_inphase`/`_outphase`,
  added the same way and in the same place as the pre-existing
  `_real`/`_imaginary`/`_ph` suffixes (same `isAddNamePostFixes` gate, no new
  CLI option).
- The JSON sidecar's `ImageType` array gained the same treatment: `WATER`/
  `FAT`/`IN_PHASE`/`OUT_OF_PHASE` are appended alongside the existing
  `MAGNITUDE`/`PHASE`/`REAL`/`IMAGINARY`/`FIELDMAPHZ` tokens, in the same
  function, using the same duplicate-guard pattern (`strstr` check before
  appending). Without this, the sidecar's `ImageType` was the generic
  dataset-level passthrough for all four Dixon outputs — identical whether
  the file was Water or Fat — so anything reading the JSON instead of the
  filename had no way to tell them apart.

## Fix 2 — stale dimension-index carryover (`// start/end dixon slice order fix`)

A pre-existing latent bug, not something introduced by Fix 1 — but Fix 1 is
what makes it visible, because it pushes the `imageType` index into the 4–7
range.

The Philips Enhanced-MR "issue 809" kludge (`isKludgeIssue809` in
`nii_dicom.cpp`) rewrites `d.dimensionIndexValues[]` once per frame:

```cpp
int d2 = d.dimensionIndexValues[2];
int d3 = d.dimensionIndexValues[3];
... d.dimensionIndexValues[3] = imageType; ...
```

`d.dimensionIndexValues[]` is a single struct-scope array reused across every
frame of the file. Each frame refreshes only indices `0 .. nDimIndxVal-1` from
the file's own `(0020,9157) DimensionIndexValues` tag. For a Dixon recon block
that only declares 3 dimensions (Stack/InStack/ImageTypeMR), `nDimIndxVal` is
3, so index 3 is **never refreshed from the file** — it silently keeps
whatever the *previous frame's* kludge run last wrote there, i.e. the
previous frame's `imageType`. That stale value (`d3`) then gets copied into
`dimIdx[6]` for the *current* frame, one frame out of phase with reality.

Downstream, `saveDcm2Nii()` decides whether to sort `dcmDim[]` forward or
reverse by comparing which array slot varies (`maxVariableItem`) against where
the file declares slice position (`stackPositionItem`). Because the stale
`d3` "varies" (in a meaningless, off-by-one-frame way), it can hijack
`maxVariableItem` away from the real `imageType` slot, picking the wrong sort
direction. `saveDcm2Nii()` then samples one representative frame per volume
with a fixed stride (`slice = i * dim3`) that assumes a clean type-major
layout — with the wrong sort direction, that stride lands mid-block, and one
volume's metadata (and therefore its Dixon type/TE/etc.) gets silently
borrowed from the wrong frame.

Fix: only read `d2`/`d3` from `d.dimensionIndexValues[]` when the *current*
frame's `nDimIndxVal` actually populated that slot; otherwise treat it as 0
(matches the file having no such dimension, rather than inheriting whatever
was last written there):

```cpp
int d2 = (nDimIndxVal > 2) ? d.dimensionIndexValues[2] : 0;
int d3 = (nDimIndxVal > 3) ? d.dimensionIndexValues[3] : 0;
```

This is gated inside the same Philips-only `isKludgeIssue809` branch as
before, so it cannot affect any other vendor, and it cannot affect Philips
files that already declare 4 dimensions (their index 3 is always freshly
written, so `d3` was never stale for them).

## Fix 3 — new fields must be zero-initialized in `clear_dicom_data()`

Any new `TDICOMdata` boolean that participates in `isSameSet()` or the
`gradDynVol` equality check (as the four Fix 1 flags do) **must** be added to
`clear_dicom_data()`'s explicit zero-initialization, alongside
`isHasPhase`/`isHasReal`/`isHasImaginary`/`isHasMagnitude`. Missing this was
caught by the regression suite: the new fields held uninitialized stack
garbage on every non-Philips file (the Philips-only detection code that would
normally set them never runs for other vendors), and since the equality
checks now compare them unconditionally, garbage values produced spurious
volume splits on unrelated vendor data. This is why the regression suite is
not optional for this kind of change — the bug was invisible on every Philips
test case and only showed up on other vendors' data.

## Fix 4 — Philips R12+ stopped populating a genuine per-echo dimension index (`// start/end echo time sort fix`)

A second, independent bug, found while investigating scrambled slice/echo
order on plain (non-Dixon) multi-echo Philips R12+ Enhanced MR data (both
gradient-echo and multi-echo spin-echo). Confirmed on paired 5.7.x/12.3.0
datasets:

```text
5.7.x  MESE: DimensionIndexValues[2] (Effective Echo Time slot) = 1,2,3,...,N   (genuine)
12.3.0 MESE: DimensionIndexValues[2] (Effective Echo Time slot) = 0,0,0,...,0   (constant!)
```

R12.3 firmware leaves the dimension index for `Effective Echo Time` at a
constant `0` for every echo — the underlying `(0018,9082) EffectiveEchoTime`
attribute itself is fine and genuinely varies, only the derived index that's
supposed to summarize it for sorting/grouping purposes has regressed.

Confirmed **not** to affect DTI/DWI: the dimension at that same array slot
for diffusion series is `Private DiffusionOrder`, which is genuinely
populated (never `0`) on both 5.7.x and 12.3.0. Diffusion b-value/gradient-
orientation numbers come from their own dedicated private tags
`(2005,1412)`/`(2005,1413)`, not the generic slot that regressed.

Consequence: inside `isKludgeIssue809`, when this slot is `0` for every
frame, nothing in the rewritten `dimIdx[]` distinguishes different echoes.
`qsort()` (used to sort `dcmDim[]`) is **not stable**, so frames tied on
every dimIdx field can come out in an unpredictable relative order — and
since `dti4D->sliceOrder[i]` (which physical frame's pixel data lands at
slot `i`) and per-index metadata (`TE[i]`, etc.) are read off the same
sorted position, the tie-break can pair the wrong frame's pixels with the
wrong slot. This is what produces a banded/interleaved appearance in
reformatted views, and can also fragment a single genuinely-coherent
multi-echo series into more output files than it should have.

Fix: when the raw echo-time slot (`d2`) is `0` but the frame's real `d.TE` is
non-zero, derive a distinguishing key from `TE` instead (microsecond
resolution, `roundf(d.TE * 1000.0f)`), so ties can no longer occur:

```cpp
if ((d2 == 0) && (!isSameFloatGE(d.TE, 0.0)))
    aslFlag = (int)roundf(d.TE * 1000.0f);
```

Scope: **not every multi-volume Philips series should be split by echo.**
Some series genuinely have no per-repeat TE anywhere in the DICOM at all
(every timing-related attribute — `EffectiveEchoTime`, the classic
`EchoTime`, Philips' own chemical-shift tags — reads as a constant `0`
across all repeats). For those, keeping the repeats as one multi-volume file
is correct; this fix only kicks in when a genuine non-zero `TE` exists to
key on, so it never fabricates a split where no real distinguishing timing
exists. It also does not gate on `bvalNum`/`gradNum` — see the code comment
at the fix site for why those can't be used to detect diffusion data here.

## Fix 5 — `isRealIsPhaseMapHz` leaking across frames within one file (`// start/end fieldmap fix`)

Found while reviewing output naming after Fix 1: a B0 field-map frame
(`Real`, legacy ImageType `B0 MAP\B0\UNSPECIFIED` — contains both "B0" and
"MAP", so the existing detection at the `ImageType` substring block is
correct) and ordinary multi-echo Real frames in the *same* Enhanced MR file
both ended up labeled `_fieldmaphz`/`_fieldmaphza`, with the true B0 map
pushed into the collision-suffixed name and the ordinary echoes claiming the
primary one.

Root cause: `d.isRealIsPhaseMapHz` is a persistent, file-scope flag that gets
set `true` once (when a B0 frame is parsed) and is **never reset per frame**
— unlike `isReal`/`isImaginary`/`isWater`/etc., which are all local variables
reset every frame alongside the `ComplexImageComponent` re-parse. So once any
frame in the file triggers it, every subsequent Real frame — B0 or not —
inherited the same flag when the per-series split assigned
`dcmList[indx].isRealIsPhaseMapHz`.

Fix: added a genuine per-frame local (`isRealIsPhaseMapHz`, matching the
persistent field it mirrors), reset every frame like its siblings, set at the
same `B0`+`MAP` detection site, and threaded
through the identical `TDCMdim` → `dcmDim[]` → `TDTI4D` → `dti4D` →
`gradDynVol` tie-check → per-series `dcmList[indx]` assignment chain used for
the Dixon flags in Fix 1. `d.isRealIsPhaseMapHz` itself is untouched, so
nothing that reads it during parsing (BIDS classification, JSON writer) is
affected — only the per-series naming/splitting decision now uses the
correct, non-leaking per-frame value.

## Isolation / blast radius

Every new code path is reached only when the data itself proves it applies:

1. `d.manufacturer == kMANUFACTURER_PHILIPS` gates all new token/flag logic.
2. The slice-order fix (Fix 2) only changes behavior inside the existing
   `isKludgeIssue809` branch (Philips Enhanced MR, software version > R10),
   and only for frames where `nDimIndxVal` is smaller than the slot being
   read — i.e. only where the pre-existing code was already reading
   undefined data.
3. The echo-time sort fix (Fix 4) only changes behavior when the raw
   dimension-index slot is exactly `0` while the frame's real `TE` is
   non-zero — verified against real DTI data that this never fires for
   genuine diffusion volumes, and it's a no-op for any file where that slot
   was never broken in the first place (pre-R12 data, or a single-echo
   series, or a series with no genuine per-repeat TE at all).
4. The fieldmap fix (Fix 5) only changes which per-series group a frame's
   `isRealIsPhaseMapHz` is read from; the underlying detection and every
   other consumer of `d.isRealIsPhaseMapHz` during parsing is untouched.
5. No new CLI options, no new `TDCMopts` fields — nothing for a user to set.
6. Non-Dixon Philips scans never set the four new booleans, so they take the
   same code path as before this change.
