import React, { useCallback, memo } from 'react';
import { Handle, Position } from 'reactflow';

interface Operation {
    value: string;
    label: string
}

export interface OperationInputParams {
    label?: string
    operations: any[]
    currentOperation?: Operation
    onChange?: (event: any) => void
    runOperation?: ({ currentOperation }: { currentOperation: Operation }) => void
    variables?: string[]
    factorOnChange?: (event: any) => void
    onChangeSetVariableValue?: (event: any) => void
    topHandle?: boolean
    bottomHandle?: boolean
}

export default memo(({ data }: { data: OperationInputParams }) => {
        const onChange = useCallback((event: any) => {
            const { onChange: dataOnChange } = data
            typeof dataOnChange === 'function' && dataOnChange(event)
        }, []);

        const changeFactor =
            useCallback((event: any) => {
                const { factorOnChange: dataFactorChange } = data
                typeof dataFactorChange === 'function' && dataFactorChange(event)
            }, []);

        let label = <label htmlFor="factorNumber" className={'mr2'}>{"Factor: "}</label>

        let input = <input id="factorNumber"
                           name="factorNumber"
                           className={'tc pv2'}
                           onChange={changeFactor}
                           type={"number"}/>


        return (
            <div className={'w5 pa3 ba br2 bw2 b--near-black bg-white-60'}>
                {data.topHandle && <Handle type="target" position={Position.Top}/>}

                <div>
                    <label htmlFor="operationSelect" className={'mr2'}>{data.label || "Operation: "}</label>
                    <select id="operationSelect"
                            name="operationSelect"
                            className={'w-100 mt2 tc pv1 nodrag'} onChange={onChange}
                            defaultValue={data.currentOperation?.value}>
                        {data.operations.map((opt: any) => <option key={opt.value} value={opt.value}>{opt.label}</option>)}
                    </select>
                </div>

                {data.factorOnChange && <div className={'mt2'}>
                        {label}
                        {input}
                    </div>}

                {data.bottomHandle && <Handle type="source" position={Position.Bottom}/>}
            </div>
        )
    }
)