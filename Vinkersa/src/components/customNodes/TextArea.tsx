import React from "react";
import {useDispatch} from "react-redux";
import {updateDataById} from "../../redux/actions";

type Props = {
    id: string,
    data: string,
    isReadOnly?: boolean
}

const TextArea: React.FC<Props> = (props: Props) => {
    const {id, data, isReadOnly} = props
    const dispatch = useDispatch()

    return (
        <textarea style={{resize: 'none', outline: 'none', border: 'none', width: '100%'}} readOnly={isReadOnly} value={data}
                  onChange={(e) => dispatch(updateDataById({id, data: e.target.value}))}/>
    )
}

export default TextArea