import { FC } from "react"
import { FileUploadButton } from "@components"

export const Toolbar: FC = () => {
  return (
    <div className="absolute top-0 z-10 flex h-20 max-h-20 w-full flex-row items-center justify-start gap-5 bg-yellow-200 px-20">
      <FileUploadButton />
    </div>
  )
}
