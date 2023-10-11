import React from "react"
import "./RunButtons.scss"
import SpeedButton from "./SpeedButton"

/**
 * Header containing buttons to control code flow execution
 *
 * @param {*} { setSpeed, runCode, reset } Function to change execution speed, function starting execution, function resetting run data
 * @returns
 */
const RunButtons = ({ setSpeed, runCode, reset }) => {
  return (
    <div className={"header"}>
      <div className="run-button" onClick={runCode}>{"Run"}</div>
      <SpeedButton setSpeed={setSpeed} speed={2000}>
        <span>{"Slow"}</span>
      </SpeedButton>
      <SpeedButton setSpeed={setSpeed} speed={1000}>
        <span>{"Medium"}</span>
      </SpeedButton>
      <SpeedButton setSpeed={setSpeed} speed={200}>
        <span>{"Fast"}</span>
      </SpeedButton>
      <div className="reset-button" onClick={reset}>{"Reset"}</div>
    </div>
  )
}

export default RunButtons