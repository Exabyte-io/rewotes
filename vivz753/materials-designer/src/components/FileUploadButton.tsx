import { ChangeEvent, FC, MouseEventHandler, useRef } from "react"

export const FileUploadButton: FC<{ handleFile: (input: string) => void }> = ({ handleFile }) => {
  const fileInput = useRef<HTMLInputElement>(null)

  const handleClick: MouseEventHandler<HTMLButtonElement> = (event) => {
    fileInput.current?.click()
    console.log()
  }

  const handleUpload = (e: ChangeEvent<HTMLInputElement>) => {
    const fileUploaded = e.target.files && e.target.files[0]

    if (fileUploaded) {
      const reader = new FileReader()
      reader.readAsText(fileUploaded, "UTF-8")
      reader.onloadend = (readerEvent: ProgressEvent<FileReader>) => {
        if (readerEvent?.target?.result) {
          handleFile(readerEvent.target.result.toString())
        }
      }
    }
  }

  return (
    <button
      className="flex items-center justify-center rounded-md bg-accent px-2 py-0.5 font-mozart text-xl uppercase tracking-widest text-light focus:outline-light"
      onClick={handleClick}
    >
      Import
      <input ref={fileInput} onChange={handleUpload} type="file" accept=".xyz,.poscar,.txt" className="hidden" />
    </button>
  )
}
