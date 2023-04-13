import { ControlGroup } from '@blueprintjs/core'
import React, { ReactNode } from 'react'

const viewHeadingStyle = {
    height: '5vh'
}

interface ViewHeadingProps {
    children: ReactNode;
}

const ViewHeading = ({children}: ViewHeadingProps) => {
  return (
    <ControlGroup style={viewHeadingStyle}>
        {children}
    </ControlGroup>

  )
}

export default ViewHeading