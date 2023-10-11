import React from "react"
import { ValueField, VariableField } from "./Inputs"

const DeclareComponent = ({ blockItem, editBlock }) => {
  return (
    <div className={"component"}>
      <VariableField value={blockItem.variable} handleChange={value => editBlock("variable", value)} />
      <span>{"="}</span>
      <ValueField value={blockItem.value} handleChange={value => editBlock("value", value)} />
    </div>
  )
}

export default DeclareComponent