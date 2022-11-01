import React, { useState, useCallback, useMemo, useEffect } from 'react';
import ReactFlow, {
    addEdge,
    FitViewOptions,
    applyNodeChanges,
    applyEdgeChanges,
    Node,
    Edge,
    NodeChange,
    EdgeChange,
    Connection,
    Background,
    Controls,
    Handle,
    useEdgesState,
    useNodesState,
} from 'reactflow';
import 'reactflow/dist/style.css';
import { Layout } from './Layout'
import NumberInput, { NumberInputParams } from './customNodes/NumberInput'
import OperationInput, { OperationInputParams } from './customNodes/OperationInput'
import TestNode, { TestNodeParams } from './customNodes/TestNode'
import AdjustmentInput, { AdjustmentInputParams } from './customNodes/AdjustmentInput'
import { BasicFlow } from "./Flows/BasicFlow";
import { makeAdjustmentInput, makeNumberInput, makeOperationInput, makeTestNode } from "./customNodes/makers";
import { Operation, operationOptions, primaryOperations, findOperation } from './lib/operations'
import { evaluateAdjustment } from "./lib/evaluateAdjustment";

interface DataVar {
    label: string
    value: string
}

const APPLICATION_TIMEOUT_MS = 5000

export const App = () => {
    /**
     *  PRIMARY APPLICATION STATE AND ASSOCIATED METHODS
     */
        // Raw Variable
    const [xVar, setXVar] = useState<DataVar>({
            label: 'x',
            value: '0',
        })
    // Wrap Raw Variable React Flow "data" prop
    const xVarInputData = {
        ...xVar,
        bottomHandle: true,
        setValue: (value: string) => setXVar((curr) => ({ ...curr, value }))
    }
    // Create XInput
    const XInput = makeNumberInput({
        id: 'xVar',
        data: xVarInputData,
        position: { x: 5, y: -35 }
    })


    const [yVar, setYVar] = useState<DataVar>({
        label: 'y',
        value: '0',
    })
    const yVarInputData: NumberInputParams = {
        ...yVar,
        topHandle: true,
        bottomHandle: true,
        setValue: (value: string) => setYVar((curr) => ({ ...curr, value }))
    }
    const YInput = makeNumberInput({
        id: 'yVar',
        data: yVarInputData,
        position: { x: 5, y: 105 },
    })


    const [primaryOperation, setPrimaryOperation] = useState<Operation>(operationOptions[2])
    const [primaryOperationOutput, setPrimaryOperationOutput] = useState("")
    const primaryOnChange = (event: any) =>
        setPrimaryOperation((curr: Operation) =>
            findOperation({ value: event.target.value }) || curr)
    const primaryOperationData = {
        label: 'Primary Operation',
        operations: primaryOperations,
        currentOperation: primaryOperation,
        topHandle: true,
        bottomHandle: true,
        onChange: primaryOnChange,
    }
    const PrimaryOp = makeOperationInput({
        id: 'primaryOperation',
        data: primaryOperationData,
        position: { x: 5, y: 255 }
    })


    const [primaryTest, setPrimaryTest] = useState<DataVar>({
        label: "Test",
        value: "output > 10"
    })
    const [primaryTestPassed, setPrimaryTestPassed] = useState("")
    const testInput = {
        ...primaryTest,
        onChange: (event: any) =>
            setPrimaryTest((curr) => ({ ...curr, value: event.target.value }))
    }
    const TestInput = makeTestNode({
        id: 'testOutput',
        data: testInput,
        position: { x: 5, y: 405 }
    })


    const [xVarAdjustOperation, setAdjustXVarOperation] = useState<Operation>(operationOptions[0])
    const xVarAdjustOperationOnChange = (event: any) =>
        setAdjustXVarOperation((curr: Operation) =>
            findOperation({ value: event.target.value }) || curr)
    const [xVarAdjustFactor, setAdjustXVarFactor] = useState("0")
    const xVarAdjust = {
        label: `adjust ${xVar.label}`,
        currentOperation: xVarAdjustOperation,
        factor: xVarAdjustFactor,
    }
    const xVarAdjustData = {
        ...xVarAdjust,
        operations: operationOptions,
        variables: [xVar, yVar],
        opOnChange: xVarAdjustOperationOnChange,
        factorOnChange: (event: any) => setAdjustXVarFactor(event.target.value),
        topHandle: 'source',
        leftHandle: 'target',
    }
    const XInputAdjust = makeAdjustmentInput({
        id: 'xVar-adjust', data: xVarAdjustData, position: { x: 365, y: 385 }
    })


    const [yVarAdjustOperation, setAdjustYVarOperation] = useState<Operation>(operationOptions[0])
    const [yVarAdjustFactor, setAdjustYVarFactor] = useState("0")
    const yVarAdjustOperationOnChange = (event: any) =>
        setAdjustYVarOperation((curr: Operation) =>
            findOperation({ value: event.target.value }) || curr)
    const yVarAdjust = {
        label: `adjust ${yVar.label}`,
        currentOperation: yVarAdjustOperation,
        factor: yVarAdjustFactor,
    }
    const yVarAdjustData = {
        ...yVarAdjust,
        operations: operationOptions,
        variables: [xVar, yVar,],
        opOnChange: yVarAdjustOperationOnChange,
        factorOnChange: (event: any) => setAdjustYVarFactor(event.target.value),
        topHandle: 'source',
        bottomHandle: 'target',
    }
    const YInputAdjust = makeAdjustmentInput({
        id: 'yVar-adjust',
        data: yVarAdjustData,
        position: { x: 365, y: 205 }
    })


    const [finalLabel, setFinalLabel] = useState('awaiting final outcome')
    const final = {
        label: finalLabel,
    }
    const FinalOutput = {
        id: 'yes-outcome',
        type: "output",
        data: final,
        position: { x: 205, y: 605 }
    }

    const initialNodes: Node[] = [
        XInput,
        YInput,
        PrimaryOp,
        TestInput,
        FinalOutput,
        XInputAdjust,
        YInputAdjust,
    ];
    const initialEdges: Edge[] = [
        { id: 'xVar-yVar', source: 'xVar', target: 'yVar' },
        { id: 'xVar-primaryOperation', source: 'yVar', target: 'primaryOperation' },
        { id: 'primaryOperation-test', source: 'primaryOperation', target: 'testOutput' },
        { id: 'test-yes', source: 'testOutput', target: 'yes-outcome' },
        { id: 'test-no', source: 'testOutput', target: 'xVar-adjust' },
        { id: 'adjust-yVar', source: 'xVar-adjust', target: 'yVar-adjust' },
        { id: 're-run', source: 'yVar-adjust', target: 'primaryOperation' },
    ];

    const [isRunning, setIsRunning] = useState(false)
    const [adjusting, setAdjusting] = useState(false)

    const [nodes, setNodes] = useNodesState<Node[]>(initialNodes);
    const [edges, setEdges] = useEdgesState<Edge[]>(initialEdges);

    // ~2 second timeout for the app anytime it is set to Run
    useEffect(() => {
        if(isRunning) {
            setFinalLabel('running 🏃')
            const timer = setTimeout(() => {
                setFinalLabel('run ended')
                setIsRunning(false)
            }, APPLICATION_TIMEOUT_MS)

            return () => clearTimeout(timer)
        }

    }, [isRunning])

    // primary test effect
    useEffect(() => {
        if(isRunning) {
            const test = primaryTest.value
            const parts = test.split(' ')

            const [targetVar, operator, suppliedVar] = parts

            const convenienceMap: { [shortcut: string]: string } = {
                output: primaryOperationOutput,
            }

            const tV = convenienceMap[targetVar]

            const allowedOperators = [
                '==',
                '!=',
                '===',
                '!==',
                '>',
                '<',
                '>=',
                '<='
            ]

            if(tV && allowedOperators.includes(operator)) {
                const test = `${tV}${operator}${+suppliedVar}`
                // the use of eval is strongly discouraged, but eval-ing this test string is fun
                const passes = eval(test)

                if(passes) {
                    setPrimaryTestPassed(passes)
                } else {
                    setAdjusting(true)
                }
            }
        }
    }, [isRunning, primaryOperationOutput])

    // PASSED!!
    useEffect(() => {
        if(primaryTestPassed) {
            const label = `test condition passed with x = ${xVar.value} and y = ${yVar.value}, the final result is ${primaryOperationOutput}`
            setIsRunning(false)
            setFinalLabel(label)
            setNodes(nds =>
                nds.map((node) => {
                    if(node.id !== 'yes-outcome') {
                        return node;
                    }
                    return {
                        ...node,
                        data: {
                            ...node.data,
                            label: label
                        },
                    };
                })
            )
        }
    }, [primaryTestPassed])

    // ADJUSTMENTS
    useEffect(() => {
        if(adjusting) {
            let operation = xVarAdjustOperation.value
            let operand = xVar.value
            let factor = xVarAdjustFactor

            const newX = evaluateAdjustment({
                operand,
                operation,
                factor,
                // setter: setXVar
            });

            setXVar((curr: DataVar) => ({ ...curr, value: newX }))

            operation = yVarAdjustOperation.value
            operand = yVar.value
            factor = yVarAdjustFactor
            const newY = evaluateAdjustment({
                operand,
                operation,
                factor,
                // setter: setYVar
            });

            setYVar((curr: DataVar) => ({ ...curr, value: newY }))

            setAdjusting(false)
        } else {
            if(isRunning) {
                // recurse...
                runPrimaryOperation()
            }
        }
    }, [adjusting])

    function runPrimaryOperation({ nodes, edges }: any = {}) {
        setIsRunning(true)

        // Do Calculation
        let output = ''
        switch (primaryOperation.value) {
            default:
                output = ''
                break

            case 'multiply':
                output = String(+xVar.value * +yVar.value)
                break

            case 'divide':
                output = String(+xVar.value / +yVar.value)
                break

            case 'add':
                output = String(+xVar.value + +yVar.value)
                break
        }

        // set primary output
        setPrimaryOperationOutput(output)
    }

    const Flow = <BasicFlow nodes={nodes} setNodes={setNodes} edges={edges} setEdges={setEdges}/>

    const calculatorState = {
        xVar,
        yVar,
        primaryOperation,
        xVarAdjust,
        yVarAdjust,
        primaryOperationOutput,
        primaryTest,
        primaryTestPassed,
        final,
    }

    const JSONRenderer = (
        <div className={'vh-75 pre pb2 bg-light-gray overflow-auto'}>
            {JSON.stringify(calculatorState, null, 2)}
        </div>
    )

    const RunButton = (
        <button className={'ba bw2 bw3-l b--white-80 ph3-ns pv3 pv2-ns bg-black-80 white-80 pointer'}
                onClick={() => runPrimaryOperation({ nodes, edges })}>
            <span className={'f3-l tracked'}>Run Operation</span>
        </button>
    )

    return <Layout JSONRenderer={JSONRenderer} FlowBuilder={Flow} RunButton={RunButton}/>
}
