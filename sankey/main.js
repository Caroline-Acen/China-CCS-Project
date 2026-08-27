// ======================================================
// main.js
// Publication-quality D3 Sankey Controller
// ======================================================

const YEARS = [2025, 2030, 2035, 2040, 2045, 2050];

const WIDTH = 700;
const HEIGHT = 460;

const LOCAL_COLOR = "#4E79A7";
const EXPORT_COLOR = "#B5B5B5";

const svgMap = {
  2025: "#svg2025",
  2030: "#svg2030",
  2035: "#svg2035",
  2040: "#svg2040",
  2045: "#svg2045",
  2050: "#svg2050",
};

//-----------------------------------------------------
// Tooltip
//-----------------------------------------------------

const tooltip = d3.select("body").append("div").attr("class", "tooltip");

//-----------------------------------------------------
// Load every year
//-----------------------------------------------------

Promise.all(YEARS.map((y) => d3.json(`data/${y}.json`))).then((allYears) => {
  allYears.forEach((data, index) => {
    drawSankey(svgMap[YEARS[index]], data, YEARS[index]);
  });
});

//-----------------------------------------------------
// Draw One Panel
//-----------------------------------------------------

function drawSankey(svgID, data, year) {
  const svg = d3.select(svgID).attr("viewBox", [0, 0, WIDTH, HEIGHT]);

  //-------------------------------------------------
  // Sankey Generator
  //-------------------------------------------------

  const sankey = d3
    .sankey()

    .nodeWidth(10)

    .nodePadding(14)

    .extent([
      [20, 20],
      [WIDTH - 20, HEIGHT - 20],
    ])

    .nodeAlign(d3.sankeyJustify)

    .iterations(64);

  //-------------------------------------------------

  const graph = sankey({
    nodes: data.nodes.map((d) => Object.assign({}, d)),

    links: data.links.map((d) => Object.assign({}, d)),
  });

  //-------------------------------------------------
  // Links
  //-------------------------------------------------

  svg
    .append("g")

    .selectAll("path")

    .data(graph.links)

    .join("path")

    .attr("class", "link")

    .attr("d", d3.sankeyLinkHorizontal())

    .attr("stroke", (d) => {
      return d.source.name === d.target.name ? LOCAL_COLOR : EXPORT_COLOR;
    })

    .attr(
      "stroke-width",

      (d) => Math.max(1, d.width),
    )

    .on("mousemove", (event, d) => {
      tooltip

        .style("opacity", 1)

        .html(
          "<b>" +
            d.source.name +
            "</b> → <b>" +
            d.target.name +
            "</b><br>" +
            d.value.toFixed(2) +
            " Gt",
        )

        .style("left", event.pageX + 15 + "px")

        .style("top", event.pageY - 20 + "px");
    })

    .on("mouseleave", () => {
      tooltip.style("opacity", 0);
    });

  //-------------------------------------------------
  // Nodes
  //-------------------------------------------------

  const node = svg
    .append("g")

    .selectAll("g")

    .data(graph.nodes)

    .join("g")

    .attr("class", "node");

  node
    .append("rect")

    .attr("x", (d) => d.x0)

    .attr("y", (d) => d.y0)

    .attr("width", (d) => d.x1 - d.x0)

    .attr(
      "height",

      (d) => Math.max(2, d.y1 - d.y0),
    )

    .attr("fill", "#444");

  //-------------------------------------------------
  // Labels
  //-------------------------------------------------

  node
    .append("text")

    .attr("x", (d) => (d.x0 < WIDTH / 2 ? d.x0 - 8 : d.x1 + 8))

    .attr("y", (d) => (d.y0 + d.y1) / 2)

    .attr("dy", "0.35em")

    .attr("text-anchor", (d) => (d.x0 < WIDTH / 2 ? "end" : "start"))

    .text((d) => d.name);

  //-------------------------------------------------
  // Optional value labels
  //-------------------------------------------------

  /*
    node.append("text")

        .attr(...)

        .text(d=>d.total.toFixed(2)+" Gt");
    */
}
