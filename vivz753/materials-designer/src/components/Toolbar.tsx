import { FileUploadButton } from "@components"
import { CrystalInput } from "@types"
import { FC, SetStateAction } from "react"

export const Toolbar: FC<{ setInput: (value: SetStateAction<CrystalInput>) => void }> = ({ setInput }) => {
  const handleFile = (input: string) => {
    setInput((prev) => {
      const newInput = { ...prev }
      newInput["crystalBasis"] = input
      return newInput
    })
  }

  return (
    <div className="absolute top-0 z-30 flex h-20 max-h-20 w-full flex-row items-center justify-start gap-5 bg-dark2 px-20">
      <FileUploadButton handleFile={handleFile} />
    </div>
  )
}
