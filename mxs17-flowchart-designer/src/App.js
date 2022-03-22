import { useCallback, useMemo, useState } from "react";

import Editor from "./components/Editor";
import Viewer from "./components/Viewer";

import Chart from "./chart-elements/Chart";

import { persistChartData, resetPersistedData } from "./helpers";
import { parseChart, evalAst } from "./parser";
import INITIAL_PANELS from "./initialState";

import "./App.css";

function App() {
  const [panels, setPanels] = useState(INITIAL_PANELS);

  const updatePanels = useCallback(
    (v) => {
      persistChartData(v.chart);
      setPanels(v);
    },
    [setPanels]
  );

  const tree = useMemo(() => {
    try {
      const [err, ast] = parseChart(panels.chart);
      if (err) return { result: err.message, astString: "{}" };
      const rawResult = evalAst(ast);
      const result =
        !rawResult && !(rawResult === false || rawResult === 0)
          ? "-"
          : rawResult;
      return { result, astString: JSON.stringify(ast, null, 2) };
    } catch (_) {
      return { result: "-", astString: "{}" };
    }
  }, [panels]);

  const resetEditor = useCallback(() => {
    resetPersistedData();
    setPanels((ps) => ({ ...ps, chart: new Chart() }));
  }, [setPanels]);

  return (
    <div className="App">
      <div className="layout">
        <Editor panels={panels} updatePanels={updatePanels} />
        <div className="json-container">
          <div className="top-panel">
            <button onClick={resetEditor} className="reset-button">
              Reset
            </button>
            {`Result: ${tree.result}`}
          </div>
          <Viewer data={tree.astString} />
        </div>
      </div>
    </div>
  );
}

export default App;
