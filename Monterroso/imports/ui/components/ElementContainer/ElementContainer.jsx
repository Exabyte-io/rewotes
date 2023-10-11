import React from "react"
import classNames from "classnames"
import { isNil } from "lodash"
import "./ElementContainer.scss"
import {v4 as uuid} from "uuid"

/**
 * ElementContainer allows elements to be dropped inside of it, styled
 *
 * @param {*} { children } list of elements to be wrapped
 */
const ElementContainer = ({ children }) => {
  return (
    <div className={classNames("drag-container")}>
      {!isNil(children) ? children.map( (child, index) => {
        return <div className={"element-item"} key={uuid()}>
          {child}
        </div>
        }) : <div></div>
      }
    </div>
  )
}

export default ElementContainer

