import React from "react"

/**
 * A controlled input element that allows standard variable names
 * underscore/letter followed by any number of letters, underscores, or numbers
 *
 * @param {*} { value, handleChange } Input value, onChange function
 * @returns
 */
const VariableField = ({ value, handleChange }) => {
  const validStarter = new RegExp(/[a-zA-Z_]/g)

  const convertToValidVariable = val => {
    for(let i = 0; i < val.length; i++) {
      if (validStarter.test(val.charAt(i))) {
        return val.slice(i)
      }
    }
    return ""
  }

  return (<input type="text" value={value} onChange={(event) => handleChange(
    convertToValidVariable(event.target.value)
  )}/>)
}

export default VariableField