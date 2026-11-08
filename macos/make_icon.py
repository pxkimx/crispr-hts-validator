"""Draw the app icon (guide-count bars, one depleted guide in red, and a check mark on a teal tile, as in
crisprhts/web/icon.svg) and write AppIcon.icns.   python make_icon.py OUT.icns   |   python make_icon.py --png OUT.png"""
import os
import subprocess
import sys
import tempfile

from PIL import Image, ImageDraw


def draw(size: int) -> Image.Image:
    s = size * 4                      # supersample, then shrink
    u = s / 64                        # the SVG is drawn on a 64-unit grid
    img = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    tile = Image.new("RGBA", (s, s))
    td = ImageDraw.Draw(tile)
    a, b = (15, 123, 116), (11, 63, 82)
    for i in range(2 * s):            # diagonal gradient, top left to bottom right
        t = i / (2 * s)
        c = tuple(int(a[k] + (b[k] - a[k]) * t) for k in range(3)) + (255,)
        td.line((i, 0, 0, i), fill=c, width=2)
    mask = Image.new("L", (s, s), 0)
    ImageDraw.Draw(mask).rounded_rectangle((3 * u, 3 * u, 61 * u, 61 * u), radius=int(13.5 * u), fill=255)
    img.paste(tile, (0, 0), mask)
    d = ImageDraw.Draw(img)

    def bar(x, top, colour):
        d.rounded_rectangle((x * u, top * u, (x + 5) * u, 52 * u), radius=int(2 * u), fill=colour)
    for x, top in [(11, 40), (18, 31), (25, 22), (32, 28), (39, 37)]:
        bar(x, top, (63, 208, 192, 255))
    bar(46, 47, (230, 103, 103, 255))
    w = int(3.4 * u)
    d.line((10 * u, 15 * u, 36 * u, 15 * u), fill=(230, 255, 251, 255), width=w)
    d.line((38 * u, 15 * u, 54 * u, 15 * u), fill=(127, 184, 255, 255), width=w)
    for x, c in [(10, (230, 255, 251)), (36, (230, 255, 251)), (38, (127, 184, 255)), (54, (127, 184, 255))]:
        d.ellipse((x * u - w / 2, 15 * u - w / 2, x * u + w / 2, 15 * u + w / 2), fill=c + (255,))
    pts = [(41 * u, 25 * u), (45 * u, 29 * u), (53 * u, 20 * u)]
    d.line(pts, fill="white", width=w, joint="curve")
    for x, y in pts:
        d.ellipse((x - w / 2, y - w / 2, x + w / 2, y + w / 2), fill="white")
    return img.resize((size, size), Image.LANCZOS)


def main(out: str) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        iconset = os.path.join(tmp, "AppIcon.iconset")
        os.mkdir(iconset)
        for base in (16, 32, 128, 256, 512):
            draw(base).save(os.path.join(iconset, f"icon_{base}x{base}.png"))
            draw(base * 2).save(os.path.join(iconset, f"icon_{base}x{base}@2x.png"))
        subprocess.run(["iconutil", "-c", "icns", iconset, "-o", out], check=True)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--png":
        draw(int(sys.argv[3]) if len(sys.argv) > 3 else 512).save(sys.argv[2])
    else:
        main(sys.argv[1] if len(sys.argv) > 1 else "AppIcon.icns")
