import { comparisonOperations, binaryOperations, unaryOperations } from "./operationDef";
// Helper function to perform an operation
const performOperation = (op, val1, val2) => {
    const operations = comparisonOperations.concat(binaryOperations).concat(unaryOperations);
    const operation = operations.filter((operation) => operation.value === op)[0];
    return operation.expression(val1, val2)
}

export default performOperation;