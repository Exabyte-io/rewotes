import React from "react"

/**
 * A controlled input element that accepts underscore, letters and numbres, and periods
 *
 * 
 * @param {*} { value, handleChange } Input value, onChange function
 * @returns
 */
const ValueField = ({ value, handleChange }) => {
  return (
    <input type="text" value={value} onChange={(event) => handleChange(
      event.target.value.replace(/\[0-9a-zA-Z_.]/g, '')
    )}/>
  )
}

export default ValueField