import React from "react"
import "./VariableComponent.scss"


/**
 * A Variable component is used to render a variable within the ComponentVisualizer, in the form of "a = 1"
 *
 * @param {*} {variable, value} The name of the variable, and it's respective value
 */
const VariableComponent = ({variable, value}) => {
  return (
    <div className={"display-container"}>
      <div className={"display-variable"}>{variable}</div>
      <span className={"display-equal"}>{"="}</span>
      <div className={"display-value"}>{value}</div>
    </div>
  )
}

export default VariableComponent