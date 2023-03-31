import React, { useCallback, memo } from 'react';
import { Handle, Position } from 'reactflow';

interface Adjustment {
    value: string;
    label: string
}

export interface AdjustmentInputParams {
    label?: string
    operations: any[]
    currentAdjustment?: Adjustment
    opOnChange?: (event: any) => void
    variables?: string[]
    factorOnChange?: (event: any) => void
    onChangeSetVariableValue?: (event: any) => void
    topHandle?: 'target' | 'source'
    leftHandle?: 'target' | 'source'
    rightHandle?: 'target' | 'source'
    bottomHandle?: 'target' | 'source'
}

export default memo(({ data }: { data: AdjustmentInputParams }) => {
    const onChange = useCallback((event: any) => {
        const { opOnChange: dataOpOnChange } = data
        typeof dataOpOnChange === 'function' && dataOpOnChange(event)
    }, []);

    const changeFactor =
        useCallback((event: any) => {
            const { factorOnChange: dataFactorChange } = data
            typeof dataFactorChange === 'function' && dataFactorChange(event)
        }, []);

    // wip -- equals
    // const onVariableSelection = useCallback((event: any) => {
    //     const newVariableValue = event.target.value
    //     const { onChangeSetVariableValue } = data
    //     onChangeSetVariableValue && onChangeSetVariableValue(newVariableValue)
    // }, [])

    let factorLabel
    let factorInput

    // wip -- equals
    // if(data.currentAdjustment?.value === 'equals') {
    //     factorLabel = <label htmlFor="variable-list" className={'mr2'}>{"Variable: "}</label>
    //     factorInput = <select id="variable-list"
    //                     name="variable-list"
    //                     className={'tc pv2'}
    //                     defaultValue={data.currentAdjustment?.value}
    //                     onChange={onVariableSelection}>
    //         {Array.isArray(data.variables) && data.variables.map((opt: any) =>
    //             <option value={opt.value} key={opt.value}>{opt.label}</option>)}
    //     </select>
    // } else {
    factorLabel = <label htmlFor="factorNumber" className={''}>{"Factor: "}</label>

    factorInput = <input id="factorNumber"
                   name="factorNumber"
                   className={'mt2 tc pv1 nodrag'}
                   onChange={changeFactor}
                   type={"number"}/>
    // }

    return (
        <div className={'w5 pa3 br2 ba b--dashed bg-white-60 shadow-1'}>
            {data.topHandle && <Handle type={data.topHandle} position={Position.Top}/>}
            {data.leftHandle && <Handle type={data.leftHandle} position={Position.Left}/>}
            {data.bottomHandle && <Handle type={data.bottomHandle} position={Position.Bottom}/>}
            {data.rightHandle && <Handle type={data.rightHandle} position={Position.Right}/>}

            <div>
                <label htmlFor="operationSelect" className={'tc mb2'}>{data.label || "Adjustment: "}</label>
                <select id="operationSelect"
                        name="operationSelect"
                        className={'w-100 tc mt2 pv1 nodrag'}
                        onChange={onChange}>
                    {data.operations.map((opt: any) => <option key={opt.value} value={opt.value}>{opt.label}</option>)}
                </select>
            </div>

            {data.factorOnChange && <div className={'mt2'}>
                {factorLabel}
                {factorInput}
            </div>
            }

        </div>
    )
})
