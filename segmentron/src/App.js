import logo from './logo.svg';
import './App.css';
import {useState} from 'react';

function Sequence(){
  const [select, setSelect] = useState(-1);
  //Segment Splits
  const splits = [700,1500,2330,3000,4000,5100,6000,6777];
  //Calculate segment indices
  let prev = 0;
  const segments = []
  const lengths = [];
  for(let i = 0; i < splits.length; i++){
    const segment = [];
    segment.push(prev+1, splits[i]);
    segments.push(segment);
    lengths.push(splits[i]-prev);
    prev = splits[i];
  }

  //const[segLengths, setSegLengths] = useState(lengths);

  const segmentBoxes = lengths.map((length,idx)=>{
    const styling = {height: `40px`, width: `${length/10}px`, background: select==idx?"#CCCCCC":"#888888", display: "inline-block", "box-sizing": "border-box", border: "1px solid #000000"};
    return <div style = {styling} onClick={() => {
      console.log(select);
      console.log(idx);
      if(select != idx){
        setSelect(idx);
      }else{
        setSelect(-1);
      }
    }}>{idx}</div>
  });
  
  return(
    <>
      <div style = {{"text-align": "center", /*"background-color": "#0666a3",*/ top: "50%", left:"50%", transform: "translate(-50%, -50%)", position: "absolute", width: "100%"}}>
        <div>
            {select<0?`Total Length: ${splits[splits.length-1]}`:`Segment ${select}, Nucleotides: (${segments[select][0]}, ${segments[select][1]})`}
        </div>
        {segmentBoxes}
      </div>
    </>
  );
}


function App() {
  return (
    <div className="App">
      <Sequence/>
    </div>
  );
}

export default App;
