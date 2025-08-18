// static/js/script.js

// JavaScript to toggle between the charts
function toggleCharts() {
    var selectedChart = document.getElementById("chartSelect").value;
    var chart1 = document.getElementById("chart1");
    var chart2 = document.getElementById("chart2");

    if (selectedChart === "chart1") {
        chart1.style.display = "block";  // Show chart1
        chart2.style.display = "none";   // Hide chart2
    } else {
        chart1.style.display = "none";   // Hide chart1
        chart2.style.display = "block";  // Show chart2
    }
}

// Additional JavaScript functions you may have used in your template
function setupRelayoutListener(graphDiv, sequences) {
    graphDiv.on('plotly_relayout', function(eventdata) {
        let start_index = Math.round(eventdata['xaxis.range[0]']);
        let end_index = Math.round(eventdata['xaxis.range[1]']);
                if (eventdata['xaxis.autorange'] != true) {
            let results = [];

            sequences.forEach(function(item) {
                let sumReads = 0;
                let sumProportion = 0;

                item.value_counts.forEach(function (entry) {
                    if (entry.index >= start_index && entry.index <= end_index) {
                        sumReads += entry.reads;
                        sumProportion += entry.proportion;
                    }
                });

                // Initialize variables for maximum value of 'reads' and its corresponding index
                let maxReads = -Infinity;
                let maxIndex = -1;

                // Iterate through the range and find the max 'reads' value
                for (let i = start_index; i < end_index; i++) {
                    if (item.value_counts[i].reads > maxReads) {
                        maxReads = item.value_counts[i].reads;
                        maxIndex = item.value_counts[i].index;
                    }
                }

                let consensus = item.value_counts[maxIndex]['consensus'];

                results.push({
                    "type": item.type,
                    "sequence": item.sequence,
                    "reads": sumReads,
                    "proportion": sumProportion.toFixed(4),
                    "consensus": consensus
                });
            });

            // Update UI with the new selection
            let title = 'Selection: ' + start_index + '-' + end_index;
            document.getElementById('header_selection').textContent = title;

            let table = document.getElementById('table_selection');
            let tbody = table.getElementsByTagName("tbody")[0];
            let tableContent = '';

            results.forEach(function(item) {
                tableContent += '<tr><td>' + item.type + '</td><td>' + item.reads + '</td><td>' + item.proportion + '</td><td class="consensus" reference=' + item.sequence + '>' + item.consensus + '</td></tr>';
            });
            tbody.innerHTML = tableContent;

            // Apply heatmap to table columns (assuming columnHeatmap is defined elsewhere)
            $('#table_selection').columnHeatmap({
                columns: [1, 2],
                inverse: true,
                colors: ['#f7f7f7', '#7fc97f', '#253494']
            });

            // Highlight matching consensus sequences
            $('.consensus').each(function() {
                let originalText = $(this).text();
                let comparedText = $(this).attr('reference');
                let highlightedText = highlightMatchingChars(originalText, comparedText);
                $(this).html(highlightedText);
            });
        } else {
                    // Handle case when autorange is reset (no selection)
                    let title = 'Selection: ' + start_index + '-' + end_index;
                    document.getElementById('header_selection').textContent = title;

                    let table = document.getElementById('table_selection');
                    let tbody = table.getElementsByTagName("tbody")[0];
                    let tableContent = '';

                    sequences.forEach(function (item) {
                        tableContent += '<tr><td>' + item.type + '</td><td>0</td><td>0</td><td></td></tr>';
                    });
                    tbody.innerHTML = tableContent;
                }
    });
}


var graph1 = {{plots['hist1'] | safe}};
var graph2 = {{plots['hist2'] | safe}};
var graphDiv1 = document.getElementById('chart1');
var graphDiv2 = document.getElementById('chart2');
var sequences = {{ sequences | safe }};
Plotly.plot('chart1',graph1,{});
Plotly.plot('chart2',graph2,{});

function highlightMatchingChars(strA, strB) {
    let result = "";
    for (let i = 0; i < strA.length; i++) {
        if (strA[i] === strB[i]) {
            result += `<span style="color: green;">${strA[i]}</span>`;
        } else {
            result += strA[i];
        }
    }
    return result;
}

$('.consensus').each(function() {
    let originalText = $(this).text();
    let comparedText = $(this).attr('reference');
    let highlightedText = highlightMatchingChars(originalText, comparedText);
    $(this).html(highlightedText);
});


    setupRelayoutListener(graphDiv1, sequences);
    setupRelayoutListener(graphDiv2, sequences);

    $('#searched_sequences').columnHeatmap({
        columns: [2, 3],
        inverse: true,
        colors: ['#f7f7f7', '#0570b0', '#d73027'] // Light gray to blue to red color scheme
    });

    $('#searched_sequences').columnHeatmap({
        columns: [4],
        inverse: false,
        colors: ['#d73027', '#f7f7f7', '#4575b4'] // Red to light gray to blue color scheme
    });

    $('.peaks-table').each(function() {
        $(this).columnHeatmap({
            columns: [3, 4],
            inverse: true,
            colors: ['#f7f7f7', '#7fc97f', '#253494'] // Light gray to green to dark blue color scheme
        });
    });

    $('#table_selection').columnHeatmap({
        columns: [1, 2],
        inverse: true,
        colors: ['#f7f7f7', '#7fc97f', '#253494'] // Light gray to green to dark blue color scheme
        // colors: ['#f7f7f7', '#fdae61', '#d73027'] // Light gray to orange to red color scheme
    });
