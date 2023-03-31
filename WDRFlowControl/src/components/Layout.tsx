import React from 'react'

interface LayoutElements {
    FlowBuilder: JSX.Element
    JSONRenderer: JSX.Element
    RunButton: JSX.Element
}

export const Layout = ({ FlowBuilder, JSONRenderer, RunButton }: LayoutElements): JSX.Element => (
    <div className={'avenir near-black'}>
        <div className={'vh-100 flex'}>
            <div className={'w-20 w-30-ns flex flex-column justify-between'}>
                <div className={'flex-grow-1 bg-light-blue pa3 br bb bw5-ns b--black-30'}>
                    <h1 className={'ma0 fw3'}>
                        Flow Chart Builder
                    </h1>

                    <p>
                        by Wesley Robinson
                    </p>

                    <div className={'mt4 tc'}>
                        {RunButton}
                    </div>
                </div>

                <div className={''}>{JSONRenderer}</div>
            </div>

            <div className={'flex-grow-1'}>{FlowBuilder}</div>
        </div>
    </div>
)