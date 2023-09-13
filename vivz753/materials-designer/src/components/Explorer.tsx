import { ChemistryIcon, TrashIcon } from "@icons"
import clsx from "clsx"
import { FC, useState } from "react"
import { v4 as uuidv4 } from "uuid"

interface Material {
  elementName: string
  elementSymbol: string
  materialId: string
}

export const Explorer: FC<{ hide: boolean; setHide: (hide: boolean) => void }> = ({ hide, setHide }) => {
  const [materialList, setMaterialList] = useState<Record<string, Material>>({})

  const deleteMaterial = (materialId: string): void => {
    // delete material from key value list
    if (materialList[materialId]) {
      setMaterialList((prevList) => {
        delete prevList[materialId]
        const newList = { ...prevList }
        return newList
      })
    }
  }

  const createMaterial = (elementName: string, elementSymbol: string): void => {
    const materialId = uuidv4()
    const material = {
      elementName,
      elementSymbol,
      materialId,
    }

    setMaterialList((prevList) => {
      prevList[materialId] = material
      const newList = { ...prevList }
      return newList
    })
  }

  return (
    <div
      className={clsx(
        hide ? "lg:w-0" : "lg:w-96",
        "smooth-transition-all absolute left-0 top-0 z-10 h-full w-1/2 pt-20 lg:relative lg:flex lg:pt-0"
      )}
    >
      <button
        onClick={() => setHide(!hide)}
        className="absolute left-0 top-0 z-20 m-2 mt-24 rounded-md bg-blue-500 px-2 py-1 lg:mt-0"
      >
        {hide ? `>>` : `<<`}
      </button>
      <div
        className={clsx(
          hide ? "-translate-x-full" : "translate-x-0",
          "smooth-transition flex h-full w-full transform flex-col gap-5 overflow-auto whitespace-nowrap rounded-md border border-white bg-blue-200 px-2"
        )}
      >
        {/* TODO: add form that takes in element name & symbol */}
        <button onClick={() => createMaterial("Lead", "Pb")} className="ml-auto mt-4 h-8 w-8 rounded-full bg-green-500">
          +
        </button>
        <p className="text-center text-xl">Explorer</p>
        <div className="flex flex-col gap-2">
          {Object.keys(materialList).map((materialId) => (
            <MaterialItem
              elementName={materialList[materialId].elementName}
              elementSymbol={materialList[materialId].elementSymbol}
              materialId={materialId}
              deleteMaterial={deleteMaterial}
            />
          ))}
        </div>
      </div>
    </div>
  )
}

const MaterialItem: FC<{
  elementName: string
  elementSymbol: string
  materialId: string
  deleteMaterial: (materialId: string) => void
}> = ({ elementName, elementSymbol, materialId, deleteMaterial }) => {
  return (
    <div className="lg:min-w-24 flex w-full flex-row items-center justify-between gap-3 rounded-md bg-black px-4 py-1.5">
      <button className="flex h-6 w-6 rounded-md">
        <ChemistryIcon className="stroke-blue-500" />
      </button>
      <div className="flex flex-col">
        <p>{elementName}</p>
        <p className="text-xs font-bold">{elementSymbol}</p>
      </div>
      <button onClick={() => deleteMaterial(materialId)} className="flex h-6 w-6 rounded-md">
        <TrashIcon className="stroke-red-500" />
      </button>
    </div>
  )
}
