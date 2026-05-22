## Randkluft – Help Guide

---

### **Quick start**
1. Upload a CSV file
2. Wait for **“File is ready for use”**
3. Click **Randkluft**
4. Review and adjust gates
5. Download results

---

## **1. Upload file**

Randkluft expects a **cell-by-marker CSV table**.

### Example input format

<table style="border-collapse: collapse; width: 100%; font-size: 0.95em;">
  <thead>
    <tr style="border-bottom: 3px solid #333;">
      <th style="text-align: right; padding: 8px 10px;"><strong>Cell</strong></th>
      <th style="text-align: right; padding: 8px 10px;"><strong>DNA1</strong></th>
      <th style="text-align: right; padding: 8px 10px;"><strong>MART1</strong></th>
      <th style="text-align: right; padding: 8px 10px;"><strong>CD207</strong></th>
      <th style="text-align: right; padding: 8px 10px;"><strong>SOX10</strong></th>
      <th style="text-align: right; padding: 8px 10px;"><strong>GZMB</strong></th>
    </tr>
  </thead>
  <tbody>
    <tr style="border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>0</strong></td>
      <td style="text-align: right; padding: 7px 10px;">17556.8360</td>
      <td style="text-align: right; padding: 7px 10px;">20322.19500</td>
      <td style="text-align: right; padding: 7px 10px;">1423.66860</td>
      <td style="text-align: right; padding: 7px 10px;">12338.48700</td>
      <td style="text-align: right; padding: 7px 10px;">2668.6628</td>
    </tr>
    <tr style="background-color: #f7f9fb; border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>1</strong></td>
      <td style="text-align: right; padding: 7px 10px;">5139.9307</td>
      <td style="text-align: right; padding: 7px 10px;">593.40510</td>
      <td style="text-align: right; padding: 7px 10px;">328.63390</td>
      <td style="text-align: right; padding: 7px 10px;">1405.91530</td>
      <td style="text-align: right; padding: 7px 10px;">794.3678</td>
    </tr>
    <tr style="border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>2</strong></td>
      <td style="text-align: right; padding: 7px 10px;">23695.3630</td>
      <td style="text-align: right; padding: 7px 10px;">396.37250</td>
      <td style="text-align: right; padding: 7px 10px;">815.00574</td>
      <td style="text-align: right; padding: 7px 10px;">581.09454</td>
      <td style="text-align: right; padding: 7px 10px;">551.8911</td>
    </tr>
    <tr style="background-color: #f7f9fb; border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>3</strong></td>
      <td style="text-align: right; padding: 7px 10px;">9178.6300</td>
      <td style="text-align: right; padding: 7px 10px;">635.32590</td>
      <td style="text-align: right; padding: 7px 10px;">357.60742</td>
      <td style="text-align: right; padding: 7px 10px;">324.48890</td>
      <td style="text-align: right; padding: 7px 10px;">798.1222</td>
    </tr>
    <tr style="border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>4</strong></td>
      <td style="text-align: right; padding: 7px 10px;">14441.8300</td>
      <td style="text-align: right; padding: 7px 10px;">455.64908</td>
      <td style="text-align: right; padding: 7px 10px;">581.14197</td>
      <td style="text-align: right; padding: 7px 10px;">822.02840</td>
      <td style="text-align: right; padding: 7px 10px;">799.6288</td>
    </tr>
    <tr style="background-color: #f7f9fb; border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>...</strong></td>
      <td style="text-align: right; padding: 7px 10px;">...</td>
      <td style="text-align: right; padding: 7px 10px;">...</td>
      <td style="text-align: right; padding: 7px 10px;">...</td>
      <td style="text-align: right; padding: 7px 10px;">...</td>
      <td style="text-align: right; padding: 7px 10px;">...</td>
    </tr>
    <tr style="border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>13195</strong></td>
      <td style="text-align: right; padding: 7px 10px;">11157.3220</td>
      <td style="text-align: right; padding: 7px 10px;">481.10248</td>
      <td style="text-align: right; padding: 7px 10px;">624.24706</td>
      <td style="text-align: right; padding: 7px 10px;">776.60986</td>
      <td style="text-align: right; padding: 7px 10px;">1207.3271</td>
    </tr>
    <tr style="background-color: #f7f9fb; border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>13196</strong></td>
      <td style="text-align: right; padding: 7px 10px;">13381.4160</td>
      <td style="text-align: right; padding: 7px 10px;">966.54550</td>
      <td style="text-align: right; padding: 7px 10px;">510.22894</td>
      <td style="text-align: right; padding: 7px 10px;">762.61110</td>
      <td style="text-align: right; padding: 7px 10px;">1138.3619</td>
    </tr>
    <tr style="border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>13197</strong></td>
      <td style="text-align: right; padding: 7px 10px;">11910.0730</td>
      <td style="text-align: right; padding: 7px 10px;">489.51560</td>
      <td style="text-align: right; padding: 7px 10px;">701.63460</td>
      <td style="text-align: right; padding: 7px 10px;">640.83850</td>
      <td style="text-align: right; padding: 7px 10px;">1175.2833</td>
    </tr>
    <tr style="background-color: #f7f9fb; border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>13198</strong></td>
      <td style="text-align: right; padding: 7px 10px;">13151.7190</td>
      <td style="text-align: right; padding: 7px 10px;">1492.03150</td>
      <td style="text-align: right; padding: 7px 10px;">795.49310</td>
      <td style="text-align: right; padding: 7px 10px;">5074.33800</td>
      <td style="text-align: right; padding: 7px 10px;">1410.3674</td>
    </tr>
    <tr style="border-bottom: 1px solid #ddd;">
      <td style="text-align: right; padding: 7px 10px;"><strong>13199</strong></td>
      <td style="text-align: right; padding: 7px 10px;">11540.0120</td>
      <td style="text-align: right; padding: 7px 10px;">436.22940</td>
      <td style="text-align: right; padding: 7px 10px;">327.97354</td>
      <td style="text-align: right; padding: 7px 10px;">2296.54420</td>
      <td style="text-align: right; padding: 7px 10px;">1149.2882</td>
    </tr>
  </tbody>
</table>

### Required columns
- **Marker intensities**
  - One or more marker columns (e.g. CD3, CD20, PD1, etc.)

### Optional columns
- **Spatial coordinates**
  - `x` and `y` centroid coordinates, if available.
  - If coordinates are not provided, Randkluft creates an internal grid for visualization.

### Important notes
- Randkluft treats each uploaded CSV as a single sample.
- DNA / Hoechst markers are **automatically masked and excluded from gating**.
- A reference input file can be downloaded from the **Upload file** tab.

---

## **2. Run Randkluft (automatic gating)**

After uploading your file:

1. Wait for the notification
   **“File is ready for use”**
   (shown as a small notice in the bottom-right corner).

2. Navigate to the **Randkluft** tab.
   - All detected markers will appear automatically.
   - By default, all markers are selected and ready for gating.

3. Review the marker list.
   - Deselect any markers you do not want to gate.

4. Click the **Randkluft** button.
   - A progress bar appears in the bottom-right corner showing
     **“Randkluft in action”** with percentage completion.

5. When gating is complete:
   - A confirmation notice appears:
     **“Gates are found!”**

---

## **3. Visual inspection and manual adjustment**

To inspect and refine gates:

- Enable **“Show gate”**.
- Marker histograms will be generated automatically.
- A notice **“Generating plot…”** appears while plots are rendered.

You can:
- Select markers to inspect
- Visually assess the estimated gates
- Manually adjust gate values using numeric input

Updated gates are applied immediately to downstream analysis.

---

## **4. Download results**

Once you are satisfied with the gates:

### Available downloads
- **Gate estimates (CSV)**
  A table containing the estimated gate value for each marker.

- **Histogram plots (PDF)**
  A PDF file containing all generated marker histograms.

- **Current plots**
  Download only the plots currently displayed in the app.

These outputs are suitable for downstream analysis, reporting, and reproducibility.

---

## **Tips**
- If gating finishes instantly without results, the marker distribution may already be negatively skewed.
- In such cases, visual inspection and manual gating are recommended.
- For best performance, ensure marker columns contain numeric values and minimal missing data.

---

For additional details, troubleshooting, and conceptual background, please consult the **FAQ** and **Contact** sections.
