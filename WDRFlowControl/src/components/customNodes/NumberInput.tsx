import React, { useCallback, memo } from 'react';
import { Handle, Position } from 'reactflow';

export interface NumberInputParams {
    label: string
    value: string
    setValue?: (value: any) => void
    topHandle?: boolean
    bottomHandle?: boolean
}

export default memo(({ data }: {data: NumberInputParams}) => {
    const onChange = useCallback((event: any) => {
        const { setValue } = data
        if (typeof setValue === 'function') setValue(event.target.value)
    }, []);

    return (
        <div className={'w5 pa3 ba br2 b--gold bg-white-60'}>
            {data.topHandle && <Handle type="target" position={Position.Top}/>}
            <div className={''}>
                <label htmlFor="number" className={''}>{data.label || "Value: "}</label>
                <input id="number"
                       name="number"
                       className={'w-100 mt2 tc pv2 nodrag'}
                       onChange={onChange}
                       type="number"/>
            </div>
            {data.bottomHandle && <Handle type="source" position={Position.Bottom}/>}
        </div>
    );
})
