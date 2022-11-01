import React, { useCallback, memo } from 'react';
import { Handle, Position } from 'reactflow';

export interface TestNodeParams {
    label?: string
    topHandle?: boolean
    bottomHandle?: boolean
    onChange?: (event: any) => void
    value?: string
}

export default memo(({ data }: { data: TestNodeParams }) => {
        const onChange = useCallback((event: any) => {
            const { onChange: dataOnChange } = data
            typeof dataOnChange === 'function' && dataOnChange(event)
        }, []);

        return (
            <div className={'w5 pa3 ba br2 b--near-black bg-white-60'}>
                {<Handle type="target" position={Position.Top}/>}
                <div>
                    <label htmlFor="function-test" className={''}>{data.label || "Value: "}</label>
                    <input id="function-test"
                           name="function-test"
                           className={'w-100 mt2 tc pv2'}
                           onChange={onChange}
                           type="text"
                           placeholder={data.value}/>
                </div>
                {<Handle type="source" position={Position.Right}/>}
            </div>
        );
    }
)