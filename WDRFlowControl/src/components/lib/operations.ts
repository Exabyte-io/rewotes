export interface Operation {
    value: string
    label: string
    factors?: string[]
    divisor?: string
    dividend?: string
    to?: string
}

export const operationOptions: Operation[] = [
    { value: "none", label: 'None' },
    { value: "increment", label: "Increment" },
    { value: "decrement", label: "Decrement" },
    { value: "multiply", label: "Multiply", factors: [], },
    { value: "divide", label: "Divide", divisor: 'x', dividend: 'y', },
    { value: "add", label: "Add" },
]

export const primaryOperations: Operation[] = [
    { value: "none", label: 'None' },
    { value: "multiply", label: "Multiply", factors: [], },
    { value: "divide", label: "Divide", divisor: 'x', dividend: 'y', },
    { value: "add", label: "Add" },
]

export const findOperation = ({ value, label }: { value?: string, label?: string }) =>
    operationOptions.find((op: Operation) => op.value === value || op.label === label)
