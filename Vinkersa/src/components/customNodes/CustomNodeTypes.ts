import { NodeTypes } from 'react-flow-renderer';
import Terminal from "./Terminal";
import IO from "./IO";

// @ts-ignore
const nodeTypes: NodeTypes = { terminal: Terminal, io: IO };

export default nodeTypes