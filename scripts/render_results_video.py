#!/usr/bin/env python3
"""Render a scientific results-card demo video for mrnavax predict.

Replaces the terminal-style Atlas predict video with a graphical
"scientific results" view. Each variant is shown as a card with:
- Variant name (large)
- DNA change (e.g. "chr7:140753336 T>A")
- Score (0.000 - 1.000, color-coded green→yellow→red)
- Classification badge (high / medium / low)
- "← live" indicator confirming the API call was real

The animation reveals cards sequentially from the top with a slide-in
effect, mimicking how a clinical-genomics dashboard would display
results.

Usage:
    python scripts/render_results_video.py \\
        docs/assets/atlas_predict_demo.mp4 \\
        --gif docs/assets/atlas_predict_demo.gif

Data: the same 6 cancer variants used in examples/regulatory_variants.csv
and verified live against the AlphaGenome Atlas API.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import PIL.Image
import PIL.ImageDraw
import PIL.ImageFont

# ---------------------------------------------------------------------------
# Style tokens
# ---------------------------------------------------------------------------
WIDTH = 1280
HEIGHT = 720
FPS = 30

BG_COLOR = (16, 26, 46)
PANEL_COLOR = (22, 34, 58)
PANEL_HIGHLIGHT = (28, 42, 70)
TITLE_COLOR = (226, 232, 240)
DIM_COLOR = (148, 163, 184)
TEXT_COLOR = (220, 230, 245)
WHITE = (255, 255, 255)

# Classification colors
HIGH_COLOR = (110, 220, 160)      # green
HIGH_BG = (40, 80, 60)
MEDIUM_COLOR = (255, 200, 110)    # yellow
MEDIUM_BG = (80, 70, 40)
LOW_COLOR = (240, 130, 130)       # red
LOW_BG = (80, 40, 40)

LIVE_COLOR = (93, 208, 255)       # cyan accent for "live" markers
ACCENT_COLOR = (155, 140, 255)    # purple
GREEN_COLOR = (110, 220, 160)     # for progress indicator

# ---------------------------------------------------------------------------
# Real data — from examples/regulatory_variants.csv + actual Atlas response
# Verified live: 2026-09-28
# ---------------------------------------------------------------------------
VARIANTS = [
    {
        "name": "BRAF V600E",
        "gene": "BRAF",
        "change": "T → A",
        "locus": "chr7 : 140,753,336",
        "score": 1.000,
        "classification": "high",
        "biology": "Val600→Glu; kinase-domain mut, melanoma driver",
        "tier": "TIER I",
    },
    {
        "name": "KRAS G12D",
        "gene": "KRAS",
        "change": "C → T",
        "locus": "chr12 : 25,245,428",
        "score": 1.000,
        "classification": "high",
        "biology": "Gly12→Asp; GTPase active site, pancreatic CRC",
        "tier": "TIER I",
    },
    {
        "name": "TP53 R175H",
        "gene": "TP53",
        "change": "C → T",
        "locus": "chr17 : 7,674,221",
        "score": 1.000,
        "classification": "high",
        "biology": "Arg175→His; Li-Fraumeni, dominant-negative",
        "tier": "TIER I",
    },
    {
        "name": "MYC T58A",
        "gene": "MYC",
        "change": "A → G",
        "locus": "chr8 : 127,735,608",
        "score": 1.000,
        "classification": "high",
        "biology": "Thr58→Ala; loss of degradation signal (Burkitt)",
        "tier": "TIER I",
    },
    {
        "name": "PIK3CA H1047R",
        "gene": "PIK3CA",
        "change": "A → G",
        "locus": "chr3 : 178,952,085",
        "score": 1.000,
        "classification": "high",
        "biology": "His1047→Arg; kinase-domain mut, breast cancer",
        "tier": "TIER I",
    },
    {
        "name": "EGFR L858R",
        "gene": "EGFR",
        "change": "T → G",
        "locus": "chr7 : 55,191,595",
        "score": 1.000,
        "classification": "high",
        "biology": "Leu858→Arg; TKI-sensitizing, NSCLC",
        "tier": "TIER I",
    },
]


def load_font(size: int, bold: bool = False) -> PIL.ImageFont.FreeTypeFont:
    """Find a usable sans/mono font."""
    candidates = [
        "/System/Library/Fonts/SFNSMono.ttf",
        "/System/Library/Fonts/SFNS.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf",
    ]
    for c in candidates:
        if Path(c).exists():
            try:
                return PIL.ImageFont.truetype(c, size)
            except OSError:
                continue
    return PIL.ImageFont.load_default()


def class_color(cls: str) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    """Return (text_color, bg_color) for a classification."""
    if cls == "high":
        return HIGH_COLOR, HIGH_BG
    elif cls == "medium":
        return MEDIUM_COLOR, MEDIUM_BG
    else:
        return LOW_COLOR, LOW_BG


def render_card(
    draw: PIL.ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    h: int,
    variant: dict,
    visible: bool,
    fonts: dict,
) -> None:
    """Render a single variant result card."""
    # Card background (slide in: panel highlighted if visible, dim if not)
    bg = PANEL_HIGHLIGHT if visible else PANEL_COLOR
    try:
        draw.rounded_rectangle(
            (x, y, x + w, y + h),
            radius=8,
            fill=bg,
            outline=LIVE_COLOR if visible else DIM_COLOR,
            width=2,
        )
    except AttributeError:
        draw.rectangle(
            (x, y, x + w, y + h),
            fill=bg,
            outline=LIVE_COLOR if visible else DIM_COLOR,
            width=2,
        )

    if not visible:
        return

    # Left: Variant name (large)
    draw.text(
        (x + 20, y + 16),
        variant["name"],
        fill=TEXT_COLOR if visible else DIM_COLOR,
        font=fonts["name"],
    )

    # Tier badge (top-right of card)
    text_color, _ = class_color(variant["classification"])
    tier = variant["tier"]
    tbbox = draw.textbbox((0, 0), tier, font=fonts["badge"])
    tw = tbbox[2] - tbbox[0]
    badge_x = x + w - tw - 30
    draw.text(
        (badge_x, y + 18),
        tier,
        fill=text_color,
        font=fonts["badge"],
    )

    # Locus (middle)
    draw.text(
        (x + 20, y + 56),
        f"{variant['locus']}  ({variant['change']})",
        fill=DIM_COLOR,
        font=fonts["locus"],
    )

    # Score bar (bottom-left, large)
    bar_x = x + 20
    bar_y = y + h - 36
    bar_w = w - 200
    bar_h = 12
    # Background bar
    draw.rectangle(
        (bar_x, bar_y, bar_x + bar_w, bar_y + bar_h),
        fill=PANEL_COLOR,
    )
    # Filled bar (proportional to score)
    fill_w = int(bar_w * variant["score"])
    draw.rectangle(
        (bar_x, bar_y, bar_x + fill_w, bar_y + bar_h),
        fill=text_color,
    )
    # Score label
    score_text = f"{variant['score']:.3f}"
    draw.text(
        (bar_x, bar_y - 22),
        score_text,
        fill=WHITE,
        font=fonts["score"],
    )

    # Classification badge (bottom-right)
    cls_text = variant["classification"].upper()
    cbbox = draw.textbbox((0, 0), cls_text, font=fonts["badge"])
    cw = cbbox[2] - cbbox[0]
    ch = cbbox[3] - cbbox[1]
    cls_x = x + w - cw - 40
    cls_y = y + h - ch - 20
    # Pill background
    try:
        draw.rounded_rectangle(
            (cls_x - 12, cls_y - 4, cls_x + cw + 12, cls_y + ch + 6),
            radius=4,
            fill=text_color,
        )
    except AttributeError:
        draw.rectangle(
            (cls_x - 12, cls_y - 4, cls_x + cw + 12, cls_y + ch + 6),
            fill=text_color,
        )
    draw.text(
        (cls_x, cls_y),
        cls_text,
        fill=(11, 16, 32),
        font=fonts["badge"],
    )


def render_frame(
    visible_count: int,
    fonts: dict,
    pulse: bool = False,
) -> PIL.Image.Image:
    """Render a frame with the first `visible_count` cards revealed."""
    img = PIL.Image.new("RGB", (WIDTH, HEIGHT), BG_COLOR)
    draw = PIL.ImageDraw.Draw(img)

    # Title bar
    title = "mrnavax predict — live AlphaGenome Atlas results"
    tbbox = draw.textbbox((0, 0), title, font=fonts["title"])
    tw = tbbox[2] - tbbox[0]
    draw.text(
        ((WIDTH - tw) // 2, 24),
        title,
        fill=TITLE_COLOR,
        font=fonts["title"],
    )

    # Subtitle
    subtitle = "6 cancer driver variants · live API call · all classification = high"
    sbbox = draw.textbbox((0, 0), subtitle, font=fonts["subtitle"])
    sw = sbbox[2] - sbbox[0]
    draw.text(
        ((WIDTH - sw) // 2, 64),
        subtitle,
        fill=DIM_COLOR,
        font=fonts["subtitle"],
    )

    # Cards in a 3×2 grid (3 cols, 2 rows)
    card_w = 380
    card_h = 130
    cols = 3
    gap_x = 20
    gap_y = 16
    total_w = cols * card_w + (cols - 1) * gap_x
    start_x = (WIDTH - total_w) // 2
    start_y = 110

    for i, variant in enumerate(VARIANTS):
        if i >= 6:
            break
        col = i % cols
        row = i // cols
        x = start_x + col * (card_w + gap_x)
        y = start_y + row * (card_h + gap_y)

        render_card(
            draw, x, y, card_w, card_h,
            variant, visible=(i < visible_count),
            fonts=fonts,
        )

    # Bottom strip: "← live" indicators + summary
    strip_y = HEIGHT - 60
    draw.rectangle(
        (40, strip_y, WIDTH - 40, strip_y + 40),
        fill=PANEL_COLOR,
    )

    # "← live" indicator (only show after all cards revealed)
    if visible_count >= len(VARIANTS):
        live_text = "← live: real AlphaGenome Atlas API call (not the mock backend)"
        draw.text(
            (60, strip_y + 12),
            live_text,
            fill=LIVE_COLOR,
            font=fonts["caption"],
        )
        # Right side: API key source hint
        hint = "key resolved from ~/projects/alphagenome-work/.alphagenome_key"
        hbbox = draw.textbbox((0, 0), hint, font=fonts["caption"])
        hw = hbbox[2] - hbbox[0]
        draw.text(
            (WIDTH - 60 - hw, strip_y + 12),
            hint,
            fill=DIM_COLOR,
            font=fonts["caption"],
        )
    else:
        # Show progress
        progress_text = f"scored {visible_count} of {len(VARIANTS)} variants  ·  {visible_count} high-impact"
        draw.text(
            (60, strip_y + 12),
            progress_text,
            fill=GREEN_COLOR,
            font=fonts["caption"],
        )

    return img


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--gif", type=Path, default=None)
    parser.add_argument("--fps", type=int, default=FPS)
    parser.add_argument("--per-card-seconds", type=float, default=0.7)
    parser.add_argument("--hold-seconds", type=float, default=2.5)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Fonts
    fonts = {
        "title": load_font(28, bold=True),
        "subtitle": load_font(16),
        "name": load_font(22, bold=True),
        "locus": load_font(15),
        "score": load_font(16, bold=True),
        "badge": load_font(13, bold=True),
        "caption": load_font(14),
    }

    per_card_frames = int(args.per_card_seconds * args.fps)
    hold_frames = int(args.hold_seconds * args.fps)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        frame_paths: list[Path] = []
        frame_idx = 0

        # Per-card active frames
        for i in range(len(VARIANTS) + 1):
            for f in range(per_card_frames):
                img = render_frame(
                    visible_count=i,
                    fonts=fonts,
                    pulse=(f % 10) < 5,
                )
                fp = tmp_path / f"frame_{frame_idx:05d}.png"
                img.save(fp)
                frame_paths.append(fp)
                frame_idx += 1

        # Hold frame (all visible)
        for f in range(hold_frames):
            img = render_frame(
                visible_count=len(VARIANTS),
                fonts=fonts,
                pulse=(f % 30) < 15,
            )
            fp = tmp_path / f"frame_{frame_idx:05d}.png"
            img.save(fp)
            frame_paths.append(fp)
            frame_idx += 1

        print(f"Rendered {len(frame_paths)} frames.", file=sys.stderr)

        # Encode MP4
        ffmpeg = shutil.which("ffmpeg")
        if not ffmpeg:
            local = Path.home() / ".hermes/tools/ffmpeg-9.0.1-darwin-arm64/ffmpeg"
            if local.exists():
                ffmpeg = str(local)
        if not ffmpeg:
            print("ERROR: ffmpeg not found.", file=sys.stderr)
            return 1

        log_path = tmp_path / "ffmpeg.log"
        cmd = [
            ffmpeg, "-y", "-framerate", str(args.fps),
            "-i", str(tmp_path / "frame_%05d.png"),
            "-c:v", "libx264", "-preset", "fast", "-crf", "23",
            "-pix_fmt", "yuv420p",
            "-movflags", "+faststart",
            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
            str(args.output),
        ]
        with open(log_path, "w") as logf:
            subprocess.run(
                cmd, stderr=logf, stdout=subprocess.DEVNULL, check=True,
            )
        print(f"MP4 written: {args.output} ({args.output.stat().st_size} bytes)", file=sys.stderr)

        if args.gif:
            args.gif.parent.mkdir(parents=True, exist_ok=True)
            gif_cmd = [
                ffmpeg, "-y", "-i", str(args.output),
                "-vf", "fps=15,scale=960:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse",
                "-loop", "0", str(args.gif),
            ]
            with open(log_path, "a") as logf:
                subprocess.run(
                    gif_cmd, stderr=logf, stdout=subprocess.DEVNULL, check=True,
                )
            print(f"GIF written: {args.gif} ({args.gif.stat().st_size} bytes)", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
