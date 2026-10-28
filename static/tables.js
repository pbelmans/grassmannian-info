/* Shared controls for the periodic tables and Freudenthal squares. */
(function () {
  "use strict";

  const fields = ["name", "dimension", "index", "euler", "type", "full", "small-message", "big-message", "fano"];
  const data = document.getElementById("atlas-isomorphisms");
  const isomorphisms = data ? JSON.parse(data.textContent) : {};

  function clearHighlight() {
    $("td.entry, circle, text").removeClass("active alternative");
  }

  $("td.entry").mouseenter(function () {
    const entry = $(this);
    const row = entry.parent();
    const label = entry.data("label");
    const dynkin = entry.data("dynkin");
    clearHighlight();

    if (dynkin) {
      const y = entry.data("y"), z = entry.data("z");
      entry.addClass("active");
      $("circle[data-label='" + dynkin + "-" + y + "'], circle[data-label='" + dynkin + "-" + z + "']", row).addClass("active");
      $("text[data-label='" + dynkin + "-" + y + "-y'], text[data-label='" + dynkin + "-" + z + "-z']", row).addClass("active");
    } else {
      const alternatives = isomorphisms[label] || [label];
      $("*[data-label='" + alternatives.join("'], *[data-label='") + "']").addClass("alternative");
      $("*[data-label='" + label + "']").removeClass("alternative").addClass("active");
    }

    fields.forEach(function (field) {
      $("td#" + field).html(entry.data(field));
    });
    $("td#full").text(entry.data("exceptional") === "known" ? "has been constructed" : "has not yet been constructed");
    const small = entry.data("small");
    $("td#small-message").text(small === "semisimple" ? "generically semisimple" :
      small === "non-semisimple" ? "generically non-semisimple" : "not known to be (non-)semisimple");
    $("td#big-message").text(entry.data("big") === "semisimple" ? "generically semisimple" : "not known to be (non-)semisimple");

    const diagram = dynkin ? $("a > svg", row) : $("svg#dynkin-" + label.split("-")[0]);
    diagram.clone().appendTo($("td#dynkin").empty());
    $("td#dynkin svg").attr("height", 35);
  });

  $("td.entry").mouseleave(function () {
    clearHighlight();
    fields.forEach(function (field) { $("td#" + field).text(""); });
    $("td#dynkin").empty();
  });

  function colors() {
    const value = $("input[name=color]:checked").val();
    if (!value) return;
    $("td.entry").removeClass("minuscule cominuscule adjoint coadjoint exceptional-known exceptional-unknown small-semisimple small-non-semisimple small-unknown big-semisimple big-unknown fano-3 fano-2 fano-1");
    $("div#typing div.form-check").removeClass("active");
    $("form#controls ul").hide();

    if (["minuscule", "cominuscule", "adjoint", "coadjoint"].includes(value)) {
      $("td.entry[data-type~=" + value + "]").addClass(value);
    }
    const options = {
      exceptional: ["exceptional-colors", ["known", "unknown"]],
      small: ["small-semisimple-colors", ["semisimple", "non-semisimple", "unknown"]],
      big: ["big-semisimple-colors", ["semisimple", "unknown"]],
      fano: ["fano-colors", ["3", "2", "1"]]
    };
    if (options[value]) {
      const [legend, states] = options[value];
      states.forEach(function (state) {
        $("td.entry[data-" + value + "=" + state + "]").addClass(value + "-" + state);
      });
      $("ul#" + legend).show();
    }
  }

  $("input[type=radio][name=color]").change(colors);
  $(function () {
    colors();
    $("input#information").change(function () {
      $("table").toggleClass("expanded", this.checked);
    });
  });
}());
