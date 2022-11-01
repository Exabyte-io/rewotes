export function evaluateAdjustment(
    {
        operand,
        operation,
        factor,
        allowedOperations = [],
        allowedFactors = [],
    }: any
): string {
    let validOperation = true
    let validFactor = true

    if(Array.isArray(allowedOperations) && allowedOperations.length > 1) {
        validOperation = allowedOperations.includes(operation)
    }

    if(Array.isArray(allowedFactors) && allowedFactors.length > 1) {
        validFactor = allowedFactors.includes(factor)
    }

    let newValue = operand
    if(validOperation && validFactor) {
        switch (operation) {
            default:
                console.log('not a valid adjustment Operation')
                break
            case "none":
                break
            case "increment":
                newValue = String(+operand + +factor)
                break
            case "decrement":
                newValue = String(+operand - +factor)
                break
            case "multiply":
                newValue = String(+operand * +factor)
                break
            case "divide":
                newValue = String(+operand / +factor)
                break
        }
    }

    return newValue
}
