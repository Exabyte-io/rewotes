import { FC, useRef, ChangeEvent, MouseEventHandler } from "react"

export const FileUploadButton: FC<{ handleFile?: (file: any) => void }> = ({ handleFile }) => {
  const fileInput = useRef<HTMLInputElement>(null)

  const handleClick: MouseEventHandler<HTMLButtonElement> = (event) => {
    fileInput.current?.click()
    console.log()
  }

  const handleUpload = (event: ChangeEvent<HTMLInputElement>) => {
    const fileUploaded = event.target.files && event.target.files[0]
    handleFile && handleFile(fileUploaded)
    // TODO: handle the file somehow
    console.log(fileUploaded)
  }

  return (
    <button
      className="flex items-center justify-center rounded-xl bg-black px-4 py-2 text-orange-500 ring ring-orange-500"
      onClick={handleClick}
    >
      Import
      {/* TODO: accept specific file type only */}
      <input ref={fileInput} onChange={handleUpload} type="file" accept="" className="hidden" />
    </button>
  )
}
