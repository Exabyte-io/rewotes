import React from "react"
import { forEach, isNil } from "lodash"
import { ElementContainer } from "../ElementContainer"
import VariableComponent from "./VariableComponent"
import { v4 as uuid } from "uuid"

/**
 *  The visualizer takes the current state of the program, and renders a human readable format within an ElementContainer component
 *
 * @param {*} { state } The current state of the program run in json format
 */
const ComponentVisualizer = ({ state }) => {
  const populateComps = (state) => {
    const comps = []
    if (!isNil(state)) {
      forEach(state.vars, (val, key) => {
        comps.push(<VariableComponent variable={key} value={val} key={uuid()}/>)
      })
      comps.push(<hr key={uuid()}/>)
      return [...populateComps(state.parentFrame), ...comps]
    }

    return comps
  }
  return (
    <ElementContainer>
      {state !== "Error" ? populateComps(state): [<div key={uuid()}>{"Error"}</div>]}
    </ElementContainer>
  )
}

export default ComponentVisualizer