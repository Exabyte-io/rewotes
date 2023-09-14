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
        hide ? "w-0" : "w-72 lg:w-96",
        "smooth-transition-all absolute left-0 top-0 z-[1] h-[calc(100vh-80px)] grow-0 text-light lg:relative lg:flex lg:pt-0"
      )}
    >
      <button
        onClick={() => setHide(!hide)}
        className="absolute bottom-0 left-0 z-[2] m-2 mt-24 whitespace-nowrap rounded-md bg-accent px-2 py-0.5 font-mozart text-xl uppercase tracking-widest text-light lg:mt-4"
      >
        {hide ? `MATERIALS >>` : `<< CLOSE`}
      </button>
      <div
        className={clsx(
          hide ? "-translate-x-full" : "translate-x-0",
          "smooth-transition flex w-full transform flex-col gap-5 overflow-auto whitespace-nowrap rounded-sm border border-dark1 bg-dark2 px-2 py-10 focus-within:border-accent"
        )}
      >
        {/* TODO: add form that takes in element name & symbol */}
        <div className="flex flex-row items-start justify-center">
          <p className="w-full px-10 text-start font-mozart text-xl uppercase tracking-widest">Materials</p>
          <button
            onClick={() => createMaterial("Lead", "Pb")}
            className="min-h-8 min-w-8 right-0 top-0 mx-4 flex shrink-0 items-center justify-center rounded-md bg-accent px-2 py-0.5 text-center font-mozart text-xl uppercase tracking-widest text-light"
          >
            <p className="flex w-full items-center justify-center text-right">{`+ ADD`}</p>
          </button>
        </div>
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
    <div className="lg:min-w-24 flex w-full flex-row items-center justify-between gap-3 rounded-md bg-dark1 px-4 py-2 text-light">
      <button className="flex h-6 w-6 rounded-md">
        <ChemistryIcon className="stroke-accent" />
      </button>
      <div className="flex w-full flex-col items-start px-4">
        <p className="font-mozart text-2xl leading-none tracking-widest">{elementName}</p>
        <p className="font-mozart text-xl font-bold leading-none tracking-widest">{elementSymbol}</p>
      </div>
      <button onClick={() => deleteMaterial(materialId)} className="flex h-6 w-6 rounded-md">
        <TrashIcon className="stroke-red-500" />
      </button>
    </div>
  )
}
