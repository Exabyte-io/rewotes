import React, {useState} from "react";

const TextArea: React.FC = () => {
    const [text, setText] = useState<string>('')

    return (
        <textarea style={{resize: 'none', outline: 'none', border: 'none', width: '100%'}} value={text}
                  onChange={(e) => setText(e.target.value)}/>
    )
}

export default TextArea