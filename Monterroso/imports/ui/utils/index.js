//Checks if the type is an Arithmetic operation
const isMathOp = val => val === "Add" || val === "Subtract" || val === "Multiply" || val === "Divide"

const checkIfNumber = val => /\d+\.?\d*$/.test(val)

export { isMathOp, checkIfNumber }