import './jsonViewer.css';

const JsonViewer = ({ jsonData }) => {
	return (
		<div id='json-viewer' className='overflow-scroll overflow-x-hidden'>
			<p>JSON Data:</p>

			<pre>{JSON.stringify(jsonData, undefined, 2)}</pre>
		</div>
	);
};

export default JsonViewer;
