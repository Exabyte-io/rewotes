import React from "react"
import { ValueField, VariableField } from "./Inputs"

const EqualityComponent = ({ blockItem, editBlock }) => {
  return (
    <div className={"component"}>
      <span>{"If"}</span>
      <VariableField value={blockItem.value1} handleChange={value => editBlock("value1", value)} />
      <span>{"==="}</span>
      <ValueField value={blockItem.value2} handleChange={value => editBlock("value2", value)} />
    </div>
  )
}

export default EqualityComponent