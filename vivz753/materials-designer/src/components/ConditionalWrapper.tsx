"use client"

import { FC, ReactElement, ReactNode } from "react"

export const ConditionalWrapper: FC<{
  condition: boolean
  wrapper: (children: ReactNode) => ReactElement<any, any> | null
  children: ReactElement<any, any> | null
}> = ({ condition, wrapper, children }) => (condition ? wrapper(children) : children)