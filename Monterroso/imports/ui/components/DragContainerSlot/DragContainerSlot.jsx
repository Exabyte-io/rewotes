import React from "react"
import classNames from "classnames"
import { v4 as uuid } from "uuid"
import "./DragContainerSlot.scss"
import { useDrop } from "react-dnd"

/**
 * A DragContainerSlot does does most of the heavy lifting of adding a dropped item into the base json.
 * It is aware of it's location within the json, and upon a drop action, inserts the added item and a new slot as required.
 * 
 * @param {*} { className, index, addToBlocks} Additional classnames, current index of the slot in the list, function to edit base json
 * @returns
 */
const DragContainerSlot = ({ className, index, addToBlocks}) => {

  const [{ isOver }, drop] = useDrop(() => ({
    accept: "Item",
    drop: (item, monitor) => {
      const slotKey = uuid()
      const itemKey = uuid()
      addToBlocks(slotKey, index + 1, {type: "Slot"})
      if (item.type === "If") {
        addToBlocks(uuid(), index + 1, {type: "Slot"})
      }
      addToBlocks(itemKey, index + 1, item)
    },
    collect: monitor => ({
      isOver: !!monitor.isOver(),
      canDrop: !!monitor.canDrop()
    })
  }))
  return (
    <div ref={drop} className={classNames("drag-container-item", {"drag-hovered": isOver })}>
    </div>
  )
}

export default DragContainerSlot

