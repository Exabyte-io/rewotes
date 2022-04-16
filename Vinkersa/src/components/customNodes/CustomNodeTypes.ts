import { NodeTypes } from 'react-flow-renderer';
import Terminal from "./Terminal";
import IO from "./IO";
import Process from "./Process";
import Decision from './Decision'

export type NodeNameTypes = 'terminal' | 'io' | 'process' | 'decision'
const nodeTypes: NodeTypes = { terminal: Terminal, io: IO, process: Process, decision: Decision };

export default nodeTypes