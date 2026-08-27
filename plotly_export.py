"""Plotly HTML export helpers (300 dpi PNG metadata on toolbar download)."""

from __future__ import annotations

import plotly.io as pio

from config import PLOTLY_PNG_DPI, PLOTLY_PNG_SCALE

# Plotly's scale upsizes pixels but the browser export omits PNG pHYs metadata,
# so Windows reports 96 dpi. Patch Plotly.downloadImage to embed pHYs at target dpi.
_DPI_DOWNLOAD_POST_SCRIPT = """
(function () {
  var TARGET_DPI = __DPI__;
  var SCALE = TARGET_DPI / 96;

  function readU32BE(bytes, offset) {
    return (
      (bytes[offset] << 24) |
      (bytes[offset + 1] << 16) |
      (bytes[offset + 2] << 8) |
      bytes[offset + 3]
    ) >>> 0;
  }

  function writeU32BE(bytes, offset, value) {
    bytes[offset] = (value >>> 24) & 0xff;
    bytes[offset + 1] = (value >>> 16) & 0xff;
    bytes[offset + 2] = (value >>> 8) & 0xff;
    bytes[offset + 3] = value & 0xff;
  }

  var crcTable = null;
  function makeCrcTable() {
    var table = new Uint32Array(256);
    for (var n = 0; n < 256; n++) {
      var c = n;
      for (var k = 0; k < 8; k++) {
        c = (c & 1) ? (0xedb88320 ^ (c >>> 1)) : (c >>> 1);
      }
      table[n] = c >>> 0;
    }
    return table;
  }

  function crc32Bytes(bytes, start, length) {
    if (!crcTable) {
      crcTable = makeCrcTable();
    }
    var c = 0xffffffff;
    for (var i = start; i < start + length; i++) {
      c = crcTable[(c ^ bytes[i]) & 0xff] ^ (c >>> 8);
    }
    return (c ^ 0xffffffff) >>> 0;
  }

  function buildPhysChunk(dpi) {
    var ppm = Math.round(dpi / 0.0254);
    var data = new Uint8Array(9);
    writeU32BE(data, 0, ppm);
    writeU32BE(data, 4, ppm);
    data[8] = 1;

    var type = new Uint8Array([0x70, 0x48, 0x59, 0x73]);
    var crcInput = new Uint8Array(13);
    crcInput.set(type, 0);
    crcInput.set(data, 4);

    var chunk = new Uint8Array(21);
    writeU32BE(chunk, 0, 9);
    chunk.set(type, 4);
    chunk.set(data, 8);
    writeU32BE(chunk, 17, crc32Bytes(crcInput, 0, 13));
    return chunk;
  }

  function embedPngDpi(arrayBuffer, dpi) {
    var bytes = new Uint8Array(arrayBuffer);
    if (bytes.length < 33) {
      return arrayBuffer;
    }

    var physChunk = buildPhysChunk(dpi);
    var parts = [bytes.slice(0, 8)];
    var pos = 8;
    var inserted = false;

    while (pos + 12 <= bytes.length) {
      var length = readU32BE(bytes, pos);
      var type = String.fromCharCode(
        bytes[pos + 4],
        bytes[pos + 5],
        bytes[pos + 6],
        bytes[pos + 7]
      );
      var end = pos + 12 + length;
      if (end > bytes.length) {
        break;
      }
      var chunk = bytes.slice(pos, end);
      if (type === "IHDR") {
        parts.push(chunk);
        parts.push(physChunk);
        inserted = true;
      } else if (type !== "pHYs") {
        parts.push(chunk);
      }
      pos = end;
    }

    if (!inserted) {
      return arrayBuffer;
    }

    var total = parts.reduce(function (sum, part) { return sum + part.length; }, 0);
    var out = new Uint8Array(total);
    var offset = 0;
    for (var i = 0; i < parts.length; i++) {
      out.set(parts[i], offset);
      offset += parts[i].length;
    }
    return out.buffer;
  }

  function downloadPngWithDpi(gd, opts) {
    opts = Object.assign({}, opts || {});
    opts.format = "png";
    opts.scale = opts.scale || SCALE;
    return Plotly.toImage(gd, opts).then(function (dataUrl) {
      return fetch(dataUrl)
        .then(function (resp) { return resp.arrayBuffer(); })
        .then(function (buffer) {
          var patched = embedPngDpi(buffer, TARGET_DPI);
          var blob = new Blob([patched], { type: "image/png" });
          var url = URL.createObjectURL(blob);
          var link = document.createElement("a");
          link.href = url;
          var filename = opts.filename || "plot";
          link.download = String(filename).replace(/[^a-z0-9_\\-\\.]/gi, "_") + ".png";
          document.body.appendChild(link);
          link.click();
          document.body.removeChild(link);
          URL.revokeObjectURL(url);
        });
    });
  }

  if (typeof Plotly === "undefined") {
    return;
  }

  Plotly.downloadImage = function (gd, opts) {
    opts = Object.assign({}, opts || {});
    if ((opts.format || "png") === "png") {
      return downloadPngWithDpi(gd, opts);
    }
    if (Plotly.toImage) {
      return Plotly.toImage(gd, opts);
    }
  };
})();
""".replace("__DPI__", str(int(PLOTLY_PNG_DPI)))


def plotly_html_config() -> dict:
    return {
        "displayModeBar": True,
        "toImageButtonOptions": {
            "format": "png",
            "scale": PLOTLY_PNG_SCALE,
        },
    }


def write_plotly_html(fig, file: str, *, config: dict | None = None) -> None:
    pio.write_html(
        fig,
        file,
        config=config or plotly_html_config(),
        include_plotlyjs="cdn",
        post_script=_DPI_DOWNLOAD_POST_SCRIPT,
    )
