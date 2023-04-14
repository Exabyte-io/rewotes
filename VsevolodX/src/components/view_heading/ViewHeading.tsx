import React, { ReactNode } from 'react'
import styles from '../view_heading/ViewHeading.module.scss';

interface ViewHeadingProps {
    children: ReactNode;
}

const ViewHeading = ({children}: ViewHeadingProps) => {
  return (
    <div className={styles.ViewHeadingStyle}>
        {children}
    </div>

  )
}

export default ViewHeading