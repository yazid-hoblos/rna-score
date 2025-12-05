import React, { useState, useEffect } from 'react';
import './App.css';

function App() {
  const [activeTab, setActiveTab] = useState('pipeline');
  const [jobId, setJobId] = useState('');
  const [jobStatus, setJobStatus] = useState('');

  return (
    <div className="App">
      <header className="header">
        <h1>🧬 RNA Score Web Interface 🧬</h1>
        <p>Extract distances, train models, and score RNA structures</p>
      </header>

      <div className="container">
        <nav className="tabs">
          <button
            className={`tab ${activeTab === 'pipeline' ? 'active' : ''}`}
            onClick={() => setActiveTab('pipeline')}
          >
            Pipeline
          </button>
          <button
            className={`tab ${activeTab === 'extract' ? 'active' : ''}`}
            onClick={() => setActiveTab('extract')}
          >
            Extract
          </button>
          <button
            className={`tab ${activeTab === 'train' ? 'active' : ''}`}
            onClick={() => setActiveTab('train')}
          >
            Train
          </button>
          <button
            className={`tab ${activeTab === 'score' ? 'active' : ''}`}
            onClick={() => setActiveTab('score')}
          >
            Score
          </button>
          <button
            className={`tab ${activeTab === 'status' ? 'active' : ''}`}
            onClick={() => setActiveTab('status')}
          >
            Status
          </button>
        </nav>

        <div className="content">
          {activeTab === 'pipeline' && <PipelineTab setJobId={setJobId} setJobStatus={setJobStatus} />}
          {activeTab === 'extract' && <ExtractTab setJobId={setJobId} setJobStatus={setJobStatus} />}
          {activeTab === 'train' && <TrainTab setJobId={setJobId} setJobStatus={setJobStatus} />}
          {activeTab === 'score' && <ScoreTab setJobId={setJobId} setJobStatus={setJobStatus} />}
          {activeTab === 'status' && <StatusTab jobId={jobId} setJobId={setJobId} />}
        </div>
      </div>
    </div>
  );
}

function ScoresTable({ scoresFile }) {
  const [scoresData, setScoresData] = useState(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    if (!scoresFile) return;

    const fetchScores = async () => {
      try {
        const response = await fetch(scoresFile.url);
        const csvText = await response.text();
        
        // Parse CSV
        const lines = csvText.trim().split('\n');
        const headers = lines[0].split(',').map(h => h.trim());
        const rows = lines.slice(1).map(line => {
          const values = line.split(',').map(v => v.trim());
          const row = {};
          headers.forEach((header, idx) => {
            row[header] = values[idx];
          });
          return row;
        });
        
        setScoresData({ headers, rows });
        setLoading(false);
      } catch (error) {
        console.error('Error fetching scores:', error);
        setLoading(false);
      }
    };

    fetchScores();
  }, [scoresFile]);

  if (loading) return <p className="loading-text">Loading scoring results...</p>;
  if (!scoresData || scoresData.rows.length === 0) {
    return <p className="loading-text">No scoring results available</p>;
  }

  return (
    <div className="scores-table-container">
      <table className="scores-table">
        <thead>
          <tr>
            {scoresData.headers.map(header => (
              <th key={header} className={header.includes('Score') ? 'score-column' : ''}>
                {header}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {scoresData.rows.map((row, idx) => (
            <tr key={idx}>
              {scoresData.headers.map(header => (
                <td key={`${idx}-${header}`} className={header.includes('Score') ? 'score-value' : ''}>
                  {header.includes('Score') ? parseFloat(row[header]).toFixed(2) : row[header]}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      <p style={{ fontSize: '0.9em', color: '#666', marginTop: '10px' }}>
        <strong>Note:</strong> Each row represents a test structure with its calculated total pseudo-energy score.
        Lower scores indicate more favorable conformations based on the training model.
      </p>
      <a href={scoresFile.url} download className="download-link" style={{ marginTop: '15px', display: 'block' }}>
        📥 Download Detailed Scores (scores.csv)
      </a>
    </div>
  );
}

function PipelineTab({ setJobId, setJobStatus }) {
  const [trainInputMode, setTrainInputMode] = useState('file'); // 'file' or 'text'
  const [scoreInputMode, setScoreInputMode] = useState('file'); // 'file' or 'text'
  const [trainFile, setTrainFile] = useState(null);
  const [scoreFile, setScoreFile] = useState(null);
  const [trainText, setTrainText] = useState('');
  const [scoreText, setScoreText] = useState('');
  const [atomMode, setAtomMode] = useState("C3'");
  const [cutoff, setCutoff] = useState(20.0);
  const [seqSep, setSeqSep] = useState(4);
  const [binWidth, setBinWidth] = useState(1.0);
  const [method, setMethod] = useState('histogram');
  const [maxScore, setMaxScore] = useState(10.0);
  const [pseudocount, setPseudocount] = useState(1e-6);
  const [includePlot, setIncludePlot] = useState(false);
  const [loading, setLoading] = useState(false);
  const [currentJobId, setCurrentJobId] = useState(null);
  const [completedJobId, setCompletedJobId] = useState(null);
  const [resultFiles, setResultFiles] = useState([]);
  const [pipelineSteps, setPipelineSteps] = useState({
    extract: 'pending',
    train: 'pending',
    score: 'pending',
    plot: 'pending'
  });

  // Poll for job status
  useEffect(() => {
    if (!loading || !currentJobId) return;

    const pollInterval = setInterval(async () => {
      try {
        const response = await fetch(`/api/job/${currentJobId}/status`);
        if (!response.ok) return;
        
        const status = await response.json();
        const message = status.message || '';

        // Update pipeline steps based on the current message from backend
        // Extract step number to determine current state
        if (message.includes('Step 1/4')) {
          setPipelineSteps({
            extract: 'running',
            train: 'pending',
            score: 'pending',
            plot: 'pending'
          });
        } else if (message.includes('Step 2/4')) {
          setPipelineSteps({
            extract: 'completed',
            train: 'running',
            score: 'pending',
            plot: 'pending'
          });
        } else if (message.includes('Step 3/4')) {
          setPipelineSteps({
            extract: 'completed',
            train: 'completed',
            score: 'running',
            plot: 'pending'
          });
        } else if (message.includes('Step 4/4')) {
          setPipelineSteps({
            extract: 'completed',
            train: 'completed',
            score: 'completed',
            plot: 'running'
          });
        } else if (message.includes('Pipeline complete')) {
          // Handle completion message
          setPipelineSteps({
            extract: 'completed',
            train: 'completed',
            score: 'completed',
            plot: 'completed'
          });
          setCompletedJobId(currentJobId);
          setLoading(false);
          clearInterval(pollInterval);
          // Fetch results
          try {
            const resultsResponse = await fetch(`/api/job/${currentJobId}/results`);
            if (resultsResponse.ok) {
              const results = await resultsResponse.json();
              setResultFiles(results.files || []);
            }
          } catch (error) {
            console.error('Error fetching results:', error);
          }
        }

        // Check if job is complete or failed (status field)
        if (status.status === 'completed') {
          setPipelineSteps({
            extract: 'completed',
            train: 'completed',
            score: 'completed',
            plot: 'completed'
          });
          setCompletedJobId(currentJobId);
          setLoading(false);
          clearInterval(pollInterval);
          // Fetch results
          try {
            const resultsResponse = await fetch(`/api/job/${currentJobId}/results`);
            if (resultsResponse.ok) {
              const results = await resultsResponse.json();
              setResultFiles(results.files || []);
            }
          } catch (error) {
            console.error('Error fetching results:', error);
          }
        } else if (status.status === 'failed') {
          const failedStep = getFailedStep(message);
          // Mark all completed steps and the failed one
          const newSteps = {
            extract: 'completed',
            train: 'completed',
            score: 'completed',
            plot: 'completed'
          };
          newSteps[failedStep] = 'failed';
          // Mark steps after failure as not started
          const stepOrder = ['extract', 'train', 'score', 'plot'];
          const failedIndex = stepOrder.indexOf(failedStep);
          for (let i = failedIndex + 1; i < stepOrder.length; i++) {
            newSteps[stepOrder[i]] = 'pending';
          }
          for (let i = 0; i < failedIndex; i++) {
            newSteps[stepOrder[i]] = 'completed';
          }
          setPipelineSteps(newSteps);
          setLoading(false);
          clearInterval(pollInterval);
        }
      } catch (error) {
        console.error('Status poll error:', error);
      }
    }, 1000); // Poll every 1 second

    return () => clearInterval(pollInterval);
  }, [loading, currentJobId]);

  const getFailedStep = (message) => {
    if (message.includes('Step 1')) return 'extract';
    if (message.includes('Step 2')) return 'train';
    if (message.includes('Step 3')) return 'score';
    if (message.includes('Step 4')) return 'plot';
    return 'extract';
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    
    // Validate training input
    let trainFileToSubmit = trainFile;
    if (trainInputMode === 'text') {
      if (!trainText.trim()) {
        alert('Please provide training PDB IDs or upload a training file');
        return;
      }
      // Create a File object from text input
      trainFileToSubmit = new File(
        [trainText],
        'training_structures.txt',
        { type: 'text/plain' }
      );
    } else if (!trainFile) {
      alert('Please select a training file');
      return;
    }

    // Validate scoring input
    let scoreFileToSubmit = scoreFile;
    if (scoreInputMode === 'text') {
      if (!scoreText.trim()) {
        alert('Please provide scoring PDB IDs or upload a scoring file');
        return;
      }
      // Create a File object from text input
      scoreFileToSubmit = new File(
        [scoreText],
        'scoring_structures.txt',
        { type: 'text/plain' }
      );
    } else if (!scoreFile) {
      alert('Please select a scoring file');
      return;
    }

    setLoading(true);
    setPipelineSteps({
      extract: 'running',
      train: 'pending',
      score: 'pending',
      plot: 'pending'
    });

    const formData = new FormData();
    formData.append('train_file', trainFileToSubmit);
    formData.append('score_file', scoreFileToSubmit);
    formData.append('atom_mode', atomMode);
    formData.append('cutoff', cutoff);
    formData.append('seq_sep', seqSep);
    formData.append('bin_width', binWidth);
    formData.append('method', method);
    formData.append('max_score', maxScore);
    formData.append('pseudocount', pseudocount);
    formData.append('include_plot', includePlot);

    try {
      console.log('Submitting pipeline request...');
      const response = await fetch('/api/pipeline', {
        method: 'POST',
        body: formData
      });
      
      console.log('Response status:', response.status, response.statusText);
      const responseText = await response.text();
      console.log('Response body:', responseText);
      
      if (!response.ok) {
        throw new Error(`HTTP ${response.status}: ${responseText}`);
      }
      
      let data;
      try {
        data = JSON.parse(responseText);
      } catch (e) {
        throw new Error(`Invalid JSON response: ${responseText}`);
      }
      
      console.log('Parsed data:', data);
      
      if (!data || !data.job_id) {
        throw new Error(`No job_id in response: ${JSON.stringify(data)}`);
      }
      
      setCurrentJobId(data.job_id);
      setJobId(data.job_id);
      setJobStatus(`Pipeline job ${data.job_id} submitted`);
      alert(`Pipeline job submitted! ID: ${data.job_id}`);
    } catch (error) {
      console.error('Pipeline submit error:', error);
      alert('Error: ' + error.message);
      setPipelineSteps(prev => ({
        ...prev,
        extract: 'failed'
      }));
      setLoading(false);
    }
  };

  return (
    <div className="pipeline-container">
      <form onSubmit={handleSubmit} className="form">
        <h2>Automated RNA Analysis Pipeline</h2>
        <p className="pipeline-description">Extract training structures → Train scoring model → Score test structures → Generate training plots</p>

        <div className="form-group">
          <label>Training Structures</label>
          <div>
            <label>
              <input
                type="radio"
                name="trainMode"
                value="file"
                checked={trainInputMode === 'file'}
                onChange={() => setTrainInputMode('file')}
                disabled={loading}
              />
              Upload File
            </label>
            <label>
              <input
                type="radio"
                name="trainMode"
                value="text"
                checked={trainInputMode === 'text'}
                onChange={() => setTrainInputMode('text')}
                disabled={loading}
              />
              Paste PDB IDs
            </label>
          </div>

          {trainInputMode === 'file' ? (
            <>
              <input
                type="file"
                onChange={(e) => setTrainFile(e.target.files[0])}
                disabled={loading}
              />
              <small>PDB/mmCIF file for extraction and training</small>
            </>
          ) : (
            <>
              <textarea
                value={trainText}
                onChange={(e) => setTrainText(e.target.value)}
                placeholder="1EHZ A&#10;1Y26 B&#10;2OEU"
                disabled={loading}
              />
              <small>Format: PDB_ID [CHAIN_ID1 CHAIN_ID2 ...], one per line</small>
            </>
          )}
        </div>

        <div className="form-group">
          <label>Test/Scoring Structures</label>
          <div>
            <label>
              <input
                type="radio"
                name="scoreMode"
                value="file"
                checked={scoreInputMode === 'file'}
                onChange={() => setScoreInputMode('file')}
                disabled={loading}
              />
              Upload File
            </label>
            <label>
              <input
                type="radio"
                name="scoreMode"
                value="text"
                checked={scoreInputMode === 'text'}
                onChange={() => setScoreInputMode('text')}
                disabled={loading}
              />
              Paste PDB IDs
            </label>
          </div>

          {scoreInputMode === 'file' ? (
            <>
              <input
                type="file"
                onChange={(e) => setScoreFile(e.target.files[0])}
                disabled={loading}
              />
              <small>PDB/mmCIF file for scoring with trained model</small>
            </>
          ) : (
            <>
              <textarea
                value={scoreText}
                onChange={(e) => setScoreText(e.target.value)}
                placeholder="1EHZ B&#10;1J8W A"
                disabled={loading}
              />
              <small>Format: PDB_ID [CHAIN_ID1 CHAIN_ID2 ...], one per line</small>
            </>
          )}
        </div>

        <fieldset className="parameter-group">
          <legend>Extraction Parameters</legend>
          
          <div className="form-row">
            <div className="form-group">
              <label>Atom Mode</label>
              <input
                type="text"
                value={atomMode}
                onChange={(e) => setAtomMode(e.target.value)}
                placeholder="e.g., C3'"
                disabled={loading}
              />
            </div>

            <div className="form-group">
              <label>Cutoff (Å)</label>
              <input
                type="number"
                value={cutoff}
                onChange={(e) => setCutoff(parseFloat(e.target.value))}
                step="0.1"
                disabled={loading}
              />
            </div>
          </div>

          <div className="form-row">
            <div className="form-group">
              <label>Sequence Separation</label>
              <input
                type="number"
                value={seqSep}
                onChange={(e) => setSeqSep(parseInt(e.target.value))}
                disabled={loading}
              />
            </div>

            <div className="form-group">
              <label>Bin Width (Å)</label>
              <input
                type="number"
                value={binWidth}
                onChange={(e) => setBinWidth(parseFloat(e.target.value))}
                step="0.1"
                disabled={loading}
              />
            </div>
          </div>

          <div className="form-group">
            <label>Method</label>
            <select value={method} onChange={(e) => setMethod(e.target.value)} disabled={loading}>
              <option value="histogram">Histogram</option>
              <option value="kde">KDE</option>
            </select>
          </div>
        </fieldset>

        <fieldset className="parameter-group">
          <legend>Training Parameters</legend>
          
          <div className="form-row">
            <div className="form-group">
              <label>Max Score</label>
              <input
                type="number"
                value={maxScore}
                onChange={(e) => setMaxScore(parseFloat(e.target.value))}
                step="0.1"
                disabled={loading}
              />
            </div>

            <div className="form-group">
              <label>Pseudocount</label>
              <input
                type="number"
                value={pseudocount}
                onChange={(e) => setPseudocount(parseFloat(e.target.value))}
                step="1e-6"
                disabled={loading}
              />
            </div>
          </div>
        </fieldset>

        <div className="form-group checkbox-group">
          <label>
            <input
              type="checkbox"
              checked={includePlot}
              onChange={(e) => setIncludePlot(e.target.checked)}
              disabled={loading}
            />
            <span>Generate training plots</span>
          </label>
        </div>

        <button type="submit" disabled={loading} className="btn-primary btn-large">
          {loading ? 'Pipeline Running...' : 'Start Pipeline'}
        </button>
      </form>

      <div className="pipeline-steps">
        <h3>Pipeline Progress</h3>
        <div className="steps-container">
          <div className={`step ${pipelineSteps.extract}`}>
            <div className="step-icon">
              {pipelineSteps.extract === 'running' && <div className="spinner"></div>}
              {pipelineSteps.extract === 'completed' && <span>✓</span>}
              {pipelineSteps.extract === 'failed' && <span>✗</span>}
              {pipelineSteps.extract === 'pending' && <span>1</span>}
            </div>
            <div className="step-label">Extract</div>
            <div className="step-status">{pipelineSteps.extract}</div>
          </div>

          <div className="step-connector"></div>

          <div className={`step ${pipelineSteps.train}`}>
            <div className="step-icon">
              {pipelineSteps.train === 'running' && <div className="spinner"></div>}
              {pipelineSteps.train === 'completed' && <span>✓</span>}
              {pipelineSteps.train === 'failed' && <span>✗</span>}
              {pipelineSteps.train === 'pending' && <span>2</span>}
            </div>
            <div className="step-label">Train</div>
            <div className="step-status">{pipelineSteps.train}</div>
          </div>

          <div className="step-connector"></div>

          <div className={`step ${pipelineSteps.score}`}>
            <div className="step-icon">
              {pipelineSteps.score === 'running' && <div className="spinner"></div>}
              {pipelineSteps.score === 'completed' && <span>✓</span>}
              {pipelineSteps.score === 'failed' && <span>✗</span>}
              {pipelineSteps.score === 'pending' && <span>3</span>}
            </div>
            <div className="step-label">Score</div>
            <div className="step-status">{pipelineSteps.score}</div>
          </div>

          <div className="step-connector"></div>

          <div className={`step ${pipelineSteps.plot}`}>
            <div className="step-icon">
              {pipelineSteps.plot === 'running' && <div className="spinner"></div>}
              {pipelineSteps.plot === 'completed' && <span>✓</span>}
              {pipelineSteps.plot === 'failed' && <span>✗</span>}
              {pipelineSteps.plot === 'pending' && <span>4</span>}
            </div>
            <div className="step-label">Plot</div>
            <div className="step-status">{pipelineSteps.plot}</div>
          </div>
        </div>
      </div>

      {completedJobId && (
        <div className="results-section">
          <h3>Pipeline Results</h3>
          {resultFiles.length > 0 ? (
            <>
              {/* Training Output Section */}
              <div className="training-results">
                <h4>Training Output</h4>
                
                {/* Training Plots */}
                {resultFiles.some(f => f.name.startsWith('plots/') && (f.name.endsWith('.png') || f.name.endsWith('.jpg'))) && (
                  <div className="results-plots">
                    <h5>Score Profiles (Training Visualizations)</h5>
                    <div className="plots-grid">
                      {resultFiles
                        .filter(f => f.name.startsWith('plots/') && (f.name.endsWith('.png') || f.name.endsWith('.jpg')))
                        .map((file) => {
                          const displayName = file.name.split('/').pop() || file.name;
                          return (
                            <div key={file.name} className="plot-container">
                              <img src={file.url} alt={displayName} className="result-plot" />
                              <p className="plot-title">{displayName}</p>
                            </div>
                          );
                        })}
                    </div>
                  </div>
                )}
                
                {/* Training Frequency/Score Tables */}
                {resultFiles.some(f => (f.name.startsWith('training_output/') || f.name.includes('freq_') || f.name.includes('score_')) && f.name.endsWith('.csv')) && (
                  <div className="training-tables">
                    <h5>Frequency & Score Tables</h5>
                    <ul className="file-list">
                      {resultFiles
                        .filter(f => (f.name.startsWith('training_output/') || f.name.includes('freq_') || f.name.includes('score_')) && f.name.endsWith('.csv'))
                        .map((file) => {
                          const displayName = file.name.split('/').pop() || file.name;
                          return (
                            <li key={file.name}>
                              <a href={file.url} download className="download-link">
                                📥 {displayName}
                              </a>
                              <span className="file-size">({(file.size / 1024).toFixed(1)} KB)</span>
                            </li>
                          );
                        })}
                    </ul>
                  </div>
                )}
              </div>

              {/* Scoring Results Section */}
              {resultFiles.some(f => f.name === 'scores.csv') && (
                <div className="scoring-results">
                  <h4>Scoring Results</h4>
                  <ScoresTable scoresFile={resultFiles.find(f => f.name === 'scores.csv')} />
                </div>
              )}

              {/* All Files Download */}
              <div className="all-files">
                <h4>All Output Files</h4>
                <ul className="file-list">
                  {resultFiles.map((file) => (
                    <li key={file.name}>
                      <a href={file.url} download className="download-link">
                        📥 {file.name}
                      </a>
                      <span className="file-size">({(file.size / 1024).toFixed(1)} KB)</span>
                    </li>
                  ))}
                </ul>
              </div>

              <div className="results-action">
                <button 
                  onClick={() => {
                    setCompletedJobId(null);
                    setResultFiles([]);
                    setPipelineSteps({
                      extract: 'pending',
                      train: 'pending',
                      score: 'pending',
                      plot: 'pending'
                    });
                  }}
                  className="btn-secondary"
                >
                  Run Another Pipeline
                </button>
              </div>
            </>
          ) : (
            <p className="loading-text">Loading results...</p>
          )}
        </div>
      )}
    </div>
  );
}

function ExtractTab({ setJobId, setJobStatus }) {
  const [file, setFile] = useState(null);
  const [atomMode, setAtomMode] = useState("C3'");
  const [distMode, setDistMode] = useState('intra');
  const [cutoff, setCutoff] = useState(20.0);
  const [seqSep, setSeqSep] = useState(4);
  const [binWidth, setBinWidth] = useState(1.0);
  const [method, setMethod] = useState('histogram');
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!file) {
      alert('Please select a file');
      return;
    }

    setLoading(true);
    const formData = new FormData();
    formData.append('file', file);
    formData.append('atom_mode', atomMode);
    formData.append('dist_mode', distMode);
    formData.append('cutoff', cutoff);
    formData.append('seq_sep', seqSep);
    formData.append('bin_width', binWidth);
    formData.append('method', method);

    try {
      const response = await fetch('/api/extract', {
        method: 'POST',
        body: formData
      });
      const data = await response.json();
      setJobId(data.job_id);
      setJobStatus(`Job ${data.job_id} submitted`);
      alert(`Job submitted! ID: ${data.job_id}`);
    } catch (error) {
      alert('Error: ' + error.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <form onSubmit={handleSubmit} className="form">
      <h2>Extract Distances</h2>
      
      <div className="form-group">
        <label>Structure File (PDB/mmCIF)</label>
        <input
          type="file"
          onChange={(e) => setFile(e.target.files[0])}
          required
        />
      </div>

      <div className="form-row">
        <div className="form-group">
          <label>Atom Mode</label>
          <input
            type="text"
            value={atomMode}
            onChange={(e) => setAtomMode(e.target.value)}
            placeholder="e.g., C3'"
          />
        </div>

        <div className="form-group">
          <label>Distance Mode</label>
          <select value={distMode} onChange={(e) => setDistMode(e.target.value)}>
            <option value="intra">Intra-chain</option>
            <option value="inter">Inter-chain</option>
          </select>
        </div>
      </div>

      <div className="form-row">
        <div className="form-group">
          <label>Cutoff (Å)</label>
          <input
            type="number"
            value={cutoff}
            onChange={(e) => setCutoff(parseFloat(e.target.value))}
            step="0.1"
          />
        </div>

        <div className="form-group">
          <label>Sequence Separation</label>
          <input
            type="number"
            value={seqSep}
            onChange={(e) => setSeqSep(parseInt(e.target.value))}
          />
        </div>

        <div className="form-group">
          <label>Bin Width (Å)</label>
          <input
            type="number"
            value={binWidth}
            onChange={(e) => setBinWidth(parseFloat(e.target.value))}
            step="0.1"
          />
        </div>
      </div>

      <div className="form-group">
        <label>Method</label>
        <select value={method} onChange={(e) => setMethod(e.target.value)}>
          <option value="histogram">Histogram</option>
          <option value="kde">KDE</option>
        </select>
      </div>

      <button type="submit" disabled={loading} className="btn-primary">
        {loading ? 'Processing...' : 'Extract'}
      </button>
    </form>
  );
}

function TrainTab({ setJobId, setJobStatus }) {
  const [files, setFiles] = useState([]);
  const [maxScore, setMaxScore] = useState(10.0);
  const [pseudocount, setPseudocount] = useState(1e-6);
  const [cutoff, setCutoff] = useState(20.0);
  const [binWidth, setBinWidth] = useState(1.0);
  const [method, setMethod] = useState('histogram');
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (files.length === 0) {
      alert('Please select histogram/KDE files');
      return;
    }

    setLoading(true);
    const formData = new FormData();
    files.forEach(file => formData.append('hist_files', file));
    formData.append('max_score', maxScore);
    formData.append('pseudocount', pseudocount);
    formData.append('cutoff', cutoff);
    formData.append('bin_width', binWidth);
    formData.append('method', method);

    try {
      const response = await fetch('/api/train', {
        method: 'POST',
        body: formData
      });
      const data = await response.json();
      setJobId(data.job_id);
      setJobStatus(`Job ${data.job_id} submitted`);
      alert(`Job submitted! ID: ${data.job_id}`);
    } catch (error) {
      alert('Error: ' + error.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <form onSubmit={handleSubmit} className="form">
      <h2>Train Scoring Function</h2>

      <div className="form-group">
        <label>Histogram/KDE Files</label>
        <input
          type="file"
          multiple
          onChange={(e) => setFiles(Array.from(e.target.files))}
          required
        />
        <small>{files.length} file(s) selected</small>
      </div>

      <div className="form-row">
        <div className="form-group">
          <label>Max Score</label>
          <input
            type="number"
            value={maxScore}
            onChange={(e) => setMaxScore(parseFloat(e.target.value))}
            step="0.1"
          />
        </div>

        <div className="form-group">
          <label>Pseudocount</label>
          <input
            type="number"
            value={pseudocount}
            onChange={(e) => setPseudocount(parseFloat(e.target.value))}
            step="1e-6"
          />
        </div>
      </div>

      <div className="form-row">
        <div className="form-group">
          <label>Cutoff (Å)</label>
          <input
            type="number"
            value={cutoff}
            onChange={(e) => setCutoff(parseFloat(e.target.value))}
            step="0.1"
          />
        </div>

        <div className="form-group">
          <label>Bin Width (Å)</label>
          <input
            type="number"
            value={binWidth}
            onChange={(e) => setBinWidth(parseFloat(e.target.value))}
            step="0.1"
          />
        </div>
      </div>

      <div className="form-group">
        <label>Method</label>
        <select value={method} onChange={(e) => setMethod(e.target.value)}>
          <option value="histogram">Histogram</option>
          <option value="kde">KDE</option>
        </select>
      </div>

      <button type="submit" disabled={loading} className="btn-primary">
        {loading ? 'Processing...' : 'Train'}
      </button>
    </form>
  );
}

function ScoreTab({ setJobId, setJobStatus }) {
  const [files, setFiles] = useState([]);
  const [cutoff, setCutoff] = useState(20.0);
  const [seqSep, setSeqSep] = useState(4);
  const [loading, setLoading] = useState(false);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (files.length === 0) {
      alert('Please select structure files');
      return;
    }

    setLoading(true);
    const formData = new FormData();
    files.forEach(file => formData.append('files', file));
    formData.append('cutoff', cutoff);
    formData.append('seq_sep', seqSep);

    try {
      const response = await fetch('/api/score', {
        method: 'POST',
        body: formData
      });
      const data = await response.json();
      setJobId(data.job_id);
      setJobStatus(`Job ${data.job_id} submitted`);
      alert(`Job submitted! ID: ${data.job_id}`);
    } catch (error) {
      alert('Error: ' + error.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <form onSubmit={handleSubmit} className="form">
      <h2>Score Structures</h2>

      <div className="form-group">
        <label>Structure Files (PDB/mmCIF)</label>
        <input
          type="file"
          multiple
          onChange={(e) => setFiles(Array.from(e.target.files))}
          required
        />
        <small>{files.length} file(s) selected</small>
      </div>

      <div className="form-row">
        <div className="form-group">
          <label>Cutoff (Å)</label>
          <input
            type="number"
            value={cutoff}
            onChange={(e) => setCutoff(parseFloat(e.target.value))}
            step="0.1"
          />
        </div>

        <div className="form-group">
          <label>Sequence Separation</label>
          <input
            type="number"
            value={seqSep}
            onChange={(e) => setSeqSep(parseInt(e.target.value))}
          />
        </div>
      </div>

      <button type="submit" disabled={loading} className="btn-primary">
        {loading ? 'Processing...' : 'Score'}
      </button>
    </form>
  );
}

function StatusTab({ jobId, setJobId }) {
  const [status, setStatus] = useState(null);
  const [results, setResults] = useState([]);
  const [loading, setLoading] = useState(false);

  const checkStatus = async () => {
    if (!jobId.trim()) {
      alert('Please enter a job ID');
      return;
    }

    setLoading(true);
    try {
      const response = await fetch(`/api/job/${jobId}/status`);
      const data = await response.json();
      setStatus(data);

      if (data.status === 'completed') {
        const resultsResponse = await fetch(`/api/job/${jobId}/results`);
        const resultsData = await resultsResponse.json();
        setResults(resultsData.files);
      }
    } catch (error) {
      alert('Error: ' + error.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="form">
      <h2>Check Job Status</h2>

      <div className="form-group">
        <label>Job ID</label>
        <input
          type="text"
          value={jobId}
          onChange={(e) => setJobId(e.target.value)}
          placeholder="Enter job ID"
        />
      </div>

      <button onClick={checkStatus} disabled={loading} className="btn-primary">
        {loading ? 'Checking...' : 'Check Status'}
      </button>

      {status && (
        <div className={`status-box status-${status.status}`}>
          <h3>Status: {status.status.toUpperCase()}</h3>
          {status.message && <p>{status.message}</p>}
          {status.error && <p className="error">{status.error}</p>}
        </div>
      )}

      {results.length > 0 && (
        <div className="results-box">
          <h3>Results</h3>
          <ul>
            {results.map((file, idx) => (
              <li key={idx}>
                <a href={file.url} download>{file.name}</a> ({(file.size / 1024).toFixed(2)} KB)
              </li>
            ))}
          </ul>
        </div>
      )}
    </div>
  );
}

export default App;
