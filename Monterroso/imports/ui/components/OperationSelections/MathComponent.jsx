import React from "react"
import { ValueField, VariableField } from "./Inputs"

/**
 * A component which stands as a mathemtical operator upon a variable
 *
 * @param {*} { blockItem, editBlock } Item in question, function to edit the item upon change
 */
const MathComponent = ({ blockItem, editBlock }) => {
  const getOperation = op => {
    if (op === "Add") {
      return "+="
    }
    else if (op === "Subtract") {
      return "-="
    }
    else if (op === "Multiply") {
      return "x="
    }
    else if (op === "Divide") {
      return "/="
    }
  }
  return (
    <div className={"component"}>
      <VariableField value={blockItem.variable} handleChange={value => editBlock("variable", value)} />
      <span>{getOperation(blockItem.type)}</span>
      <ValueField value={blockItem.value} handleChange={value => editBlock("value", value)} />
    </div>
  )
}

export default MathComponent