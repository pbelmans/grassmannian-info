// Entry point for the CodeMirror 6 browser bundle.
// Everything is imported from a single npm source so there are no
// duplicate-instance problems between @codemirror/* peer packages.
export { basicSetup, EditorView } from "codemirror";
export { EditorState } from "@codemirror/state";
export { keymap } from "@codemirror/view";
export { javascript } from "@codemirror/lang-javascript";
