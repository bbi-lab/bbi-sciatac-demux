/********** Define components **********/

/*
** General functions.
*/
function Header(props) {
  return React.createElement(
    "nav",
    { className: "navbar navbar-expand-md sticky-top navbar-light", style: { backgroundColor: "#e3f2fd" } },
    React.createElement(
      "div",
      { className: "navbar-collapse collapse w-100 order-1 order-md-0 dual-collapse2" },
      React.createElement(
        "ul",
        { className: "navbar-nav mr-auto" },
        React.createElement("img", { src: "img/bbi_icon.png", height: "70", className: "d-inline-block align-top", alt: "" })
      )
    ),
    React.createElement(
      "div",
      { className: "mx-auto order-0" },
      React.createElement(
        "a",
        { className: "navbar-brand mx-auto", href: "#" },
        "Demultiplexing ",
        props.run_name,
        " QC Dashboard"
      )
    ),
    React.createElement("div", { className: "navbar-collapse collapse w-100 order-3 dual-collapse2" })
  );
}

function TitleRow(props) {
  return React.createElement(
    "th",
    { scope: "col" },
    props.samp
  );
}

function RegRow(props) {
  return React.createElement(
    "td",
    null,
    props.val
  );
}

function NamedRow(props) {
  return React.createElement(
    "tr",
    null,
    React.createElement(
      "th",
      { scope: "row" },
      props.name
    ),
    React.createElement(
      "td",
      null,
      props.val
    )
  );
}

/*
** Sample barcode functions.
*/
function SampleTab(props) {
  var bad_wells_exist = Object.entries(props.bad_wells_barcodes).length === 0 && props.bad_wells_barcodes.constructor === Object;
  var plates = plate_list.map(function (plate_name, index) {
    return React.createElement(SampleTabPlate, { key: index, lane: props.lane, plate_name: plate_name, norm: props.norm });
  });
  return React.createElement(
    "div",
    { className: props.className, id: "navsample-lane" + props.lane, role: "tabpanel", "aria-labelledby": "navsample-lane" + props.lane + "-tab" },
    plates,
    bad_wells_exist ? React.createElement("span", null) : React.createElement(BadWellTable, { bad_wells: props.bad_wells, lane: props.lane, bad_wells_barcodes: props.bad_wells_barcodes })
  );
}

function SampleTabPlate(props) {
  var norm = props.norm;
  return React.createElement(
    "span",
    null,
    React.createElement(
      "h4",
      null,
      "Plate ",
      props.plate_name,
      " file: ",
      "img/L00" + props.lane + "_" + props.plate_name + ".tag_i5_plate.png"
    ),
    React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".tag_i5_plate.png", width: "50%", className: "rounded mx-auto d-block", alt: "..." }),
    norm ? React.createElement(
      "span",
      null,
      React.createElement(
        "h5",
        null,
        "Sentinel Normalized file: ",
        "img/L00" + props.lane + "_" + props.plate_name + ".tag_i5_plate_sent_norm.png"
      ),
      React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".tag_i5_plate_sent_norm.png", width: "50%", className: "rounded mx-auto d-block", alt: "..." }),
      React.createElement(
        "h5",
        null,
        "Barnyard Normalized file: ",
        "img/L00" + props.lane + "_" + props.plate_name + ".tag_i5_plate_barn_norm.png"
      ),
      React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".tag_i5_plate_barn_norm.png", width: "50%", className: "rounded mx-auto d-block", alt: "..." })
    ) : React.createElement("span", null)
  );
}

/* Make 'Sample Barcodes' header above lane tabs. */
function SampleBarcodes(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "sample", role: "tabpanel", "aria-labelledby": "sample-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "sample" },
        "Sample Barcodes"
      )
    ),
    React.createElement(
      "nav",
      null,
      React.createElement(
        "div",
        { className: "nav nav-tabs", id: "navsample-tab", role: "tablist" },
        props.sample_tab_head
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-sampleContent" },
      props.sample_tabs
    )
  );
}

function BadWellTable(props) {
  var Lane = "Lane " + props.lane;
  return React.createElement(
    "span",
    null,
    React.createElement(
      "h4",
      null,
      "Wells outside plate layout for " + Lane
    ),
    React.createElement(
      "table",
      { className: "table table-hover" },
      React.createElement(
        "thead",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement("th", { scope: "col" }),
          React.createElement(
            "th",
            { scope: "col" },
            "ReadCount"
          )
        )
      ),
      React.createElement(
        "tbody",
        null,
        props.bad_wells_barcodes.map(function (item, index) {
          return React.createElement(NamedRow, { key: index, name: item, val: props.bad_wells[[Lane]][[item]] });
        })
      )
    )
  );
}

/*
** PCR barcode functions.
*/
function PCRTab(props) {
  var plates = pcr_combo_list.map(function (plate_name, index) {
    return React.createElement(PCRTabPlate, { key: index, lane: props.lane, plate_name: plate_name });
  });
  return React.createElement(
    "div",
    { className: props.className, id: "navpcr-lane" + props.lane, role: "tabpanel",
      "aria-labelledby": "navpcr-lane" + props.lane + "-tab" },
    React.createElement(
      "p",
      null,
      "Wells circled in ",
      React.createElement(
        "strong",
        null,
        React.createElement(
          "span",
          { style: { color: "red" } },
          "red"
        )
      ),
      " have >2 orders of magnitude fewer reads than the median well in the plate."
    ),
    plates
  );
}

function PCRTabPlate(props) {
  return React.createElement(
    "span",
    null,
    React.createElement(
      "h4",
      null,
      "PCR Combo ",
      props.plate_name
    ),
    React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".pcr_plate.png", width: "50%",
      className: "rounded mx-auto d-block", alt: "..." })
  );
}

/* Make 'PCR Barcodes' header above lane tabs. */
function PCRBarcodes(props) {
  var lane_list = props.lane_stats.map(function (item) {
    return item.Lane;
  });
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "pcr", role: "tabpanel", "aria-labelledby": "pcr-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "pcr" },
        "PCR Barcodes"
      )
    ),
    Object.keys(props.pcr_combo_list).length === 0 ? React.createElement(
      "span",
      null,
      React.createElement(
        "h4",
        null,
        "Well combination read counts:"
      ),
      React.createElement(
        "table",
        { className: "table table-hover" },
        React.createElement(
          "thead",
          null,
          React.createElement(
            "tr",
            null,
            React.createElement(
              "th",
              { scope: "col" },
              "p5 well"
            ),
            React.createElement(
              "th",
              { scope: "col" },
              "p7 well"
            ),
            lane_list.map(function (item, index) {
              return React.createElement(TitleRow, { key: index, samp: item });
            })
          )
        ),
        React.createElement(
          "tbody",
          null,
          props.pcr_well_info.map(function (item, index) {
            return React.createElement(
              "tr",
              { key: index },
              React.createElement(
                "th",
                { scope: "row" },
                item.p5
              ),
              React.createElement(
                "th",
                { scope: "row" },
                item.p7
              ),
              lane_list.map(function (lane, index) {
                return React.createElement(RegRow, { key: index, val: item[lane] });
              })
            );
          })
        )
      )
    ) : React.createElement(
      "span",
      null,
      React.createElement(
        "nav",
        null,
        React.createElement(
          "div",
          { className: "nav nav-tabs", id: "navpcr-tab", role: "tablist" },
          props.pcr_tab_head
        )
      ),
      React.createElement(
        "div",
        { className: "tab-content", id: "nav-pcrContent" },
        props.pcr_tabs
      )
    )
  );
}

/*
** Ligation barcode functions.
*/
function LigTab(props) {
  var plates = lig_combo_list.map(function (plate_name, index) {
    return React.createElement(LigTabPlate, { key: index, lane: props.lane, plate_name: plate_name });
  });
  return React.createElement(
    "div",
    { className: props.className, id: "navlig-lane" + props.lane, role: "tabpanel", "aria-labelledby": "navlig-lane" + props.lane + "-tab" },
    plates
  );
}

function LigTabPlate(props) {
  return React.createElement(
    "span",
    null,
    React.createElement(
      "h4",
      null,
      "Ligation Plate ",
      props.plate_name
    ),
    React.createElement("img", { src: "img/L00" + props.lane + "_" + props.plate_name + ".tag_i7_plate.png", width: "50%", className: "rounded mx-auto d-block", alt: "..." })
  );
}

/* Make 'Ligation Barcodes header above lane tabs. */
function LigBarcodes(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "lig", role: "tabpanel", "aria-labelledby": "lig-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "lig-name" },
        "Ligation Barcodes"
      )
    ),
    React.createElement(
      "nav",
      null,
      React.createElement(
        "div",
        { className: "nav nav-tabs", id: "navlig-tab", role: "tablist" },
        props.lig_tab_head
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-ligContent" },
      props.lig_tabs
    )
  );
}

/*
** Recovery page functions.
*/
function RecoveryLog(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade", id: "log", role: "tabpanel", "aria-labelledby": "log-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3 border-bottom" },
      React.createElement(
        "h1",
        { className: "h3", id: "log-name" },
        "Recovery Summary"
      )
    ),
    React.createElement(
      "nav",
      null,
      React.createElement(
        "div",
        { className: "nav nav-tabs", id: "navlig-tab", role: "tablist" },
        props.log_tab_head
      )
    ),
    React.createElement(
      "div",
      { className: "tab-content", id: "nav-ligContent" },
      props.log_tabs
    )
  );
}

/*
** Log page functions.
*/
function CodeChunk(props) {
  return React.createElement(
    "pre",
    { style: { paddingLeft: '20px' } },
    React.createElement(
      "code",
      null,
      '\n' + props.text + '\n\n'
    )
  );
}

function LogTab(props) {
  return React.createElement(
    "div",
    { className: props.className, id: "navlog-lane" + props.lane, role: "tabpanel", "aria-labelledby": "navlog-lane" + props.lane + "-tab" },
    React.createElement(CodeChunk, { text: props.log })
  );
}

/*
** Summary table functions.
*/
function SummaryTable(props) {
  return React.createElement(
    "div",
    { className: "tab-pane fade show active", id: "summary", role: "tabpanel", "aria-labelledby": "summary-tab" },
    React.createElement(
      "div",
      { className: "d-flex justify-content-between flex-wrap flex-md-nowrap align-items-center pt-3 pb-2 mb-3" },
      React.createElement(
        "h1",
        { className: "h3", id: "summary-name" },
        "Summary Table"
      )
    ),
    React.createElement(
      "table",
      { className: "table table-hover" },
      React.createElement(
        "thead",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement("th", { scope: "col" }),
          props.lane_stats.map(function (item, index) {
            return React.createElement(TitleRow, { key: index, samp: item.Lane });
          })
        )
      ),
      React.createElement(
        "tbody",
        null,
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            "Total input reads"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: parseFloat(item.total_input_reads) });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            "Total corrected barcode sets"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: parseFloat(item.all_barcodes) });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            "Total corrected barcode sets not in samplesheet"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: parseFloat(item.total_not_specified_in_samplesheet) });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            "Total passed"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: parseFloat(item.all_barcodes) - parseFloat(item.total_not_specified_in_samplesheet) });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row" },
            "Pass percentage"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: (100.0 * (parseFloat(item.all_barcodes) - parseFloat(item.total_not_specified_in_samplesheet)) / parseFloat(item.total_input_reads)).toFixed(2) });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row", title: "Percent of reads where one of the barcode pieces was uncorrectable" },
            "Percent bad barcode sets"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: (100.0 * (parseFloat(item.total_input_reads) - parseFloat(item.all_barcodes)) / parseFloat(item.total_input_reads)).toFixed(2) });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row", title: "Percent of reads where the tagmentation barcodes were not correctable" },
            "Percent tagmentation barcodes uncorrectable"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: (100.0 * (parseFloat(item.total_input_reads) - parseFloat(item.tagmentation)) / parseFloat(item.total_input_reads)).toFixed(2) });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row", title: "Percent of reads where the PCR barcodes were not correctable" },
            "Percent PCR barcodes uncorrectable"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: (100.0 * (parseFloat(item.total_input_reads) - parseFloat(item.pcr)) / parseFloat(item.total_input_reads)).toFixed(2) });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row", title: "Percent of reads where the tagmentation barcodes mismatch" },
            "Percent tagmentation barcode pair mismatch"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: (100.0 * (parseFloat(item.tagmentation) - parseFloat(item.tagmentation_match)) / parseFloat(item.total_input_reads)).toFixed(2) });
          })
        ),
        React.createElement(
          "tr",
          null,
          React.createElement(
            "th",
            { scope: "row", title: "Percent of reads where the PCR barcodes mismatch" },
            "Percent PCR barcode pair mismatch"
          ),
          props.lane_stats.map(function (item, index) {
            return React.createElement(RegRow, { key: index, val: (100.0 * (parseFloat(item.pcr) - parseFloat(item.pcr_match)) / parseFloat(item.total_input_reads)).toFixed(2) });
          })
        )
      )
    )
  );
}

/*
** Demux page functions.
*/
function DemuxPage(props) {
  return React.createElement(
    "span",
    null,
    React.createElement(Header, { run_name: props.run_name }),
    React.createElement(
      "div",
      { className: "container-fluid" },
      React.createElement(
        "div",
        { className: "row" },
        React.createElement(
          "nav",
          { className: "col-md-2 d-none d-md-block bg-light sidebar" },
          React.createElement(
            "div",
            { className: "sidebar-sticky" },
            React.createElement(
              "div",
              { className: "nav flex-column nav-pills", id: "v-pills-tab", role: "tablist", "aria-orientation": "vertical" },
              React.createElement(
                "a",
                { className: "nav-link active", id: "summary-tab", "data-toggle": "pill", href: "#summary", role: "tab", "aria-controls": "summary", "aria-selected": "true" },
                "Summary Table"
              ),
              React.createElement(
                "a",
                { className: "nav-link", id: "sample-tab", "data-toggle": "pill", href: "#sample", role: "tab", "aria-controls": "sample", "aria-selected": "false" },
                "Sample Barcodes"
              ),
              React.createElement(
                "a",
                { className: "nav-link", id: "pcr-tab", "data-toggle": "pill", href: "#pcr", role: "tab", "aria-controls": "pcr", "aria-selected": "false" },
                "PCR Barcodes"
              ),
              props.level == 3 ? React.createElement(
                "a",
                { className: "nav-link", id: "lig-tab", "data-toggle": "pill", href: "#lig", role: "tab", "aria-controls": "lig", "aria-selected": "false" },
                "Ligation Barcodes"
              ) : '',
              React.createElement(
                "a",
                { className: "nav-link", id: "log-tab", "data-toggle": "pill", href: "#log", role: "tab", "aria-controls": "log", "aria-selected": "false" },
                "Recovery Summary"
              )
            )
          )
        ),
        React.createElement(
          "main",
          { role: "main", className: "col-md-9 ml-sm-auto col-lg-10 px-4", style: { paddingTop: "15px" } },
          React.createElement(
            "div",
            { className: "tab-content", id: "nav-tabContent" },
            React.createElement(SummaryTable, { lane_stats: props.lane_stats }),
            React.createElement(SampleBarcodes, { sample_tabs: props.sample_tabs, bad_wells: props.bad_wells, sample_tab_head: props.sample_tab_head }),
            React.createElement(PCRBarcodes, { pcr_tabs: props.pcr_tabs, lane_stats: props.lane_stats,
              pcr_combo_list: props.pcr_combo_list, pcr_well_info: props.pcr_well_info,
              pcr_tab_head: props.pcr_tab_head }),
            props.level == 3 ? React.createElement(LigBarcodes, { lig_tabs: props.lig_tabs, lig_tab_head: props.lig_tab_head }) : "",
            React.createElement(RecoveryLog, { log_tabs: props.log_tabs, log_tab_head: props.log_tab_head })
          )
        )
      )
    )
  );
}

/********** Generate pages **********/

/*
** run_data values are read from run_data.js
*/
var lane_list = run_data['lane_list'];
var lane_stats = run_data['lane_stats'];
var run_name = run_data['run_name'];
var plate_list = run_data['plate_list'];
var pcr_combo_list = run_data['pcr_combo_list'];
var lig_combo_list = run_data['lig_combo_list'];
var level = run_data['level'];
var bad_wells = run_data['bad_wells'];
var bad_wells_barcodes = run_data['bad_wells_barcodes'];
var include_norm = run_data['include_norm'];
var pcr_well_info = run_data['pcr_well_info'];

var LigTabs = lane_list.map(function (lane, index) {
  return React.createElement(LigTab, { key: index, className: "tab-pane fade", lane: lane });
});

var LogTabs = lane_list.map(function (lane, index) {
  return React.createElement(LogTab, { key: index, className: "tab-pane fade", lane: lane, log: log_data[lane] });
});

var SampleTabs = lane_list.map(function (lane, index) {
  return React.createElement(SampleTab, { key: index, className: "tab-pane fade", lane: lane, bad_wells: bad_wells, bad_wells_barcodes: bad_wells_barcodes, norm: include_norm });
});

var PCRLaneTabs = lane_list.map(function (lane, index) {
  return Object.keys(pcr_combo_list).length === 0 ? "" : React.createElement(PCRTab, { key: index, className: "tab-pane fade", lane: lane });
});

var LigTabHead = lane_list.map(function (lane, index) {
  return React.createElement(
    "a",
    { key: index, className: "nav-item nav-link", id: "navlig-lane" + lane + "-tab", "data-toggle": "tab", href: "#navlig-lane" + lane, role: "tab", "aria-controls": "navlig-lane" + lane, "aria-selected": "false" },
    "Lane " + lane
  );
});

var SampleTabHead = lane_list.map(function (lane, index) {
  return React.createElement(
    "a",
    { key: index, className: "nav-item nav-link", id: "navsample-lane" + lane + "-tab", "data-toggle": "tab", href: "#navsample-lane" + lane, role: "tab", "aria-controls": "navsample-lane" + lane, "aria-selected": "false" },
    "Lane " + lane
  );
});

var PCRTabHead = lane_list.map(function (lane, index) {
  return React.createElement(
    "a",
    { key: index, className: "nav-item nav-link", id: "navpcr-lane" + lane + "-tab", "data-toggle": "tab", href: "#navpcr-lane" + lane, role: "tab", "aria-controls": "navpcr-lane" + lane, "aria-selected": "false" },
    "Lane " + lane
  );
});

var LogTabHead = lane_list.map(function (lane, index) {
  return React.createElement(
    "a",
    { key: index, className: "nav-item nav-link", id: "navlog-lane" + lane + "-tab", "data-toggle": "tab", href: "#navlog-lane" + lane, role: "tab", "aria-controls": "navlog-lane" + lane, "aria-selected": "false" },
    "Lane " + lane
  );
});

ReactDOM.render(React.createElement(DemuxPage, { sample_tabs: SampleTabs, pcr_tabs: PCRLaneTabs, pcr_well_info: pcr_well_info,
  pcr_combo_list: pcr_combo_list, lig_tabs: LigTabs, lane_stats: lane_stats,
  level: level, run_name: run_name, bad_wells: bad_wells, sample_tab_head: SampleTabHead,
  lig_tab_head: LigTabHead, pcr_tab_head: PCRTabHead, log_tabs: LogTabs, log_tab_head: LogTabHead }), document.getElementById('demux_page'));