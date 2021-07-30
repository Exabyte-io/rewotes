import classNames from "classnames"
import React from "react" 
import { useDrag } from "react-dnd"
import "./DragComponent.scss"


/**
 * A DragComponent is allows itself to be dragged and be inserted into the DragContainer, where the specified item will be added to the state
 *
 * @param {*} { item, className, children } The item to be added to state, additional className, and children to be rendered within the component
 * @returns
 */
const DragComponent = ({ item, className, children }) => {
  const [, drag] = useDrag(() => ({
    type: "Item",
    item
  }))
  return (
    <div ref={drag} className={ classNames("drag-component", className)}>
      {children}
    </div>
  )
}

export default DragComponent