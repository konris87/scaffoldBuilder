#include "ScaffoldGenerator.h"
#include "Visualize/visualize.h"
#include "Wall/Wall.h"
#include "Optimization/objective.h"
#include "Optimization/bfgs.h"
#include <vtkAppendPolyData.h>
#include <vtkNamedColors.h>
#include <vtkContourFilter.h>
#include <vtkSignedDistance.h>
#include <vtkPolyDataNormals.h>
#include <vtkImplicitModeller.h>
#include <vtkCleanPolyData.h>
#include <vtkSTLWriter.h>
#include <vtkKdTreePointLocator.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkmContour.h>
#include <voro++.hh>
#include <Eigen/Dense>
#include <Utils/Utils.h>
#include <chrono>
#include <omp.h>
#include <limits>

float ScaffoldGenerator::estimate_connectivity() {

	connectedVertices.clear();
	connectedIndices.clear();

	// use the graph to estimate the percentage of connected pores (voronoi seeds) compared to domain volume
	if (!graph) {
		std::cerr << "Graph is not initialized. Please run populate_voro() first." << std::endl;
		return -1.0f;
	}

	// get a deep copy of the graph
	auto graphCopy = std::make_unique<Graph>(*graph);

	//float limitEdge = thickness;

	//graphCopy->remove_edges_above(pullbackRatio);
	graphCopy->remove_edges_below(thickness * 0.5);
	//graphCopy->remove_edges_below(1);

	// get number of vertices after removal
	int vertexNr = graphCopy->get_vertex_count();

	if (vertexNr == 0) {
		std::cerr << " Graph has zero vertices " << std::endl;
		return -1.0f;
	}

	int maxLength = 0;
	maxLength = graphCopy->find_longest_network();
	float connectivity_percentage =
		(static_cast<float>(maxLength) / vertexNr) * 100.0f;

	connectedVertices.reserve(seeds.size() * 3);
	for (const auto& s : seeds) {
		connectedVertices.push_back(s[0]);
		connectedVertices.push_back(s[1]);
		connectedVertices.push_back(s[2]);
	}

	// the centroid ids
	int centroidId = static_cast<int>(seeds.size());

	// map (seedA,seedB) -> centroid vertex index
	std::unordered_map<long long, int> pairToCentroidId;
	auto pack = [](int a, int b)->long long {
		if (a > b) std::swap(a, b);
		return (static_cast<long long>(a) << 32) | static_cast<unsigned int>(b);
	};

	const auto& adjList = graphCopy->get_adj_list();

	for (const auto& [v1, edges] : adjList) {
		for (const auto& [v2, w] : edges) {
			if (v1 >= v2) continue; // undirected, handle once

			// find (v1,v2) centroid
			Eigen::Vector3d c;
			bool found = false;

			// centroids is std::unordered_map<int, std::unordered_map<int, Eigen::Vector3d>>
			// You stored only with seedId < neighborId. We’re iterating with v1 < v2 (due to continue above),
			// so we should look up centroids[v1][v2].
			auto it1 = centroids.find(v1);
			if (it1 != centroids.end()) {
				auto it2 = it1->second.find(v2);
				if (it2 != it1->second.end()) {
					c = it2->second;
					found = true;
				}
			}

			if (!found) {
				// Fallback: midpoint of seeds (still ok for viz)
				Eigen::Vector3d p1(seeds[v1][0], seeds[v1][1], seeds[v1][2]);
				Eigen::Vector3d p2(seeds[v2][0], seeds[v2][1], seeds[v2][2]);
				c = 0.5 * (p1 + p2);
			}

			// get or create centroid vertex index
			long long key = pack(v1, v2);
			auto itc = pairToCentroidId.find(key);
			int cId;

			if (itc == pairToCentroidId.end()) {
				// append centroid position to buffer
				connectedVertices.push_back(c.x());
				connectedVertices.push_back(c.y());
				connectedVertices.push_back(c.z());

				cId = centroidId++;
				pairToCentroidId.emplace(key, cId);
			}
			else {
				cId = itc->second;
			}

			// edges: v1 -- cId -- v2
			connectedIndices.push_back(static_cast<unsigned int>(v1));
			connectedIndices.push_back(static_cast<unsigned int>(cId));
			connectedIndices.push_back(static_cast<unsigned int>(cId));
			connectedIndices.push_back(static_cast<unsigned int>(v2));
		}
	}

	return connectivity_percentage;
};

void ScaffoldGenerator::generate_mesh(
	const double& thickness,
	vtkSmartPointer<vtkPolyData>& finalPolyData,
	const std::array<int, 3>& res
) {

	// first get the center of the scaffold mesh
	double center[3] = { 0.0, 0.0, 0.0 };

	vtkIdType numPoints = scaffoldMesh->GetNumberOfPoints();

	vtkSmartPointer<vtkPoints> points = scaffoldMesh->GetPoints();
	for (vtkIdType i = 0; i < numPoints; i++) {
		double p[3];
		points->GetPoint(i, p);
		center[0] += p[0];
		center[1] += p[1];
		center[2] += p[2];
	}
	center[0] /= numPoints;
	center[1] /= numPoints;
	center[2] /= numPoints;

	std::cout << "found center" << std::endl;

	// Transform the polydata
	vtkNew<vtkTransform> transform;
	transform->Translate(center[0], center[1], center[2]);
	float sf = 1.0f -  (thickness + 0.1)  * 0.1f;
	sf = 1.0f;
	transform->Scale(sf, sf, sf);
	transform->Translate(-center[0], -center[1], -center[2]);

	vtkNew<vtkTransformPolyDataFilter> transformPD;
	transformPD->SetTransform(transform);
	transformPD->SetInputData(scaffoldMesh);
	transformPD->Update();
	scaffoldMesh->DeepCopy(transformPD->GetOutput());

	vtkNew<vtkNamedColors> colors;

	// decide mesh resolution
	int dim[3] = {};

	if (res.empty()) {
		if (thickness < 0.3) {
			//dim[0] = 100;
			//dim[1] = 100;
			//dim[2] = 100;
			dim[0] = 50;
			dim[1] = 50;
			dim[2] = 50;
		}
		else if (0.3 <= thickness && thickness < 0.5) {
			dim[0] = 75;
			dim[1] = 75;
			dim[2] = 75;
		}
		else {
			dim[0] = 100;
			dim[1] = 100;
			dim[2] = 100;
		}
	}
	else {
		// If resolution was provided, use it
		dim[0] = res[0];
		dim[1] = res[1];
		dim[2] = res[2];
	}

	//std::cout << "Res: " << dim[0] << " " << dim[1] << " " << dim[2] << std::endl;
	
	double boundsLocal[6];
	scaffoldMesh->GetBounds(boundsLocal);
	double dx = boundsLocal[1] - boundsLocal[0];
	double dy = boundsLocal[3] - boundsLocal[2];
	double dz = boundsLocal[5] - boundsLocal[4];
	double maxDim = std::max({ dx, dy, dz });

	double minFeature = thickness / 4.0; // conservative default
	// optional: compute minimal edge length across mesh (expensive) and replace minFeature

	// spacing: want a voxel spacing much smaller than minFeature
	double spacing = std::max(minFeature * 0.25, thickness * 0.125); // tuneable
	int recommendedDim = static_cast<int>(std::ceil(maxDim / spacing)) + 1;
	recommendedDim = std::clamp(recommendedDim, 32, 400); // clamps to sane ranges
	
	//std::cout << "Res: " << recommendedDim << " " << recommendedDim << " " << recommendedDim << std::endl;

	std::cout << "entering implict" << std::endl;
	auto start_time = std::chrono::steady_clock::now();

	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->AddInputData(scaffoldMesh);
	implictModeller->SetMaximumDistance(thickness / 2 * 0.9);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	//implictModeller->SetSampleDimensions(recommendedDim, recommendedDim, recommendedDim);
	implictModeller->SetModelBounds(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
	implictModeller->Update();

	auto end_time = std::chrono::steady_clock::now();

	// Calculate the duration in milliseconds
	auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	// Format the output string
	std::ostringstream oss;
	std::cout << "implicit modeler finished in "
		<< std::fixed << std::setprecision(3) // Set precision to 3 decimal places
		<< duration_ms.count() / 1000.0   // Convert ms to seconds
		<< " seconds." << std::endl;
	
	start_time = std::chrono::steady_clock::now();

	//extract isosurface using the GPU filter
	vtkNew<vtkmContour> isoFilter;
	isoFilter->SetInputConnection(implictModeller->GetOutputPort());
	isoFilter->SetValue(0, thickness / 2.0f);
	isoFilter->Update();

	end_time = std::chrono::steady_clock::now();

	// Calculate the duration in milliseconds
	duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	// Format the output string
	std::cout << "contour m filter finished in "
		<< std::fixed << std::setprecision(3) // Set precision to 3 decimal places
		<< duration_ms.count() / 1000.0   // Convert ms to seconds
		<< " seconds." << std::endl;

	// update normals
	vtkNew<vtkPolyDataNormals> norms;
	norms->SetInputConnection(isoFilter->GetOutputPort());
	norms->ComputePointNormalsOn();
	norms->ComputeCellNormalsOff();
	norms->ConsistencyOn();
	norms->AutoOrientNormalsOn();
	//norms->Update();

	// clean duplicated points
	vtkNew<vtkCleanPolyData> cleaner;
	cleaner->SetInputConnection(norms->GetOutputPort());
	cleaner->Update();

	// force triangulation (removes strips, polys with >3 vertices, etc.)
	vtkNew<vtkTriangleFilter> triFilter;
	triFilter->SetInputConnection(cleaner->GetOutputPort());
	triFilter->Update();

	// final clean to drop verts/lines
	vtkNew<vtkCleanPolyData> finalClean;
	finalClean->SetInputConnection(triFilter->GetOutputPort());
	finalClean->Update();

	vtkSmartPointer<vtkPolyData> poly = finalClean->GetOutput();

	// explicitly clear verts/lines just in case
	poly->SetVerts(nullptr);
	poly->SetLines(nullptr);
	poly->SetStrips(nullptr);

	// assign to output
	finalPolyData = poly;

	//// Debug: check composition
	//std::cout << "Final polydata stats:\n"
	//	<< "  Points: " << finalPolyData->GetNumberOfPoints() << "\n"
	//	<< "  Cells: " << finalPolyData->GetNumberOfCells() << "\n"
	//	<< "  Polys: " << finalPolyData->GetNumberOfPolys() << "\n"
	//	<< "  Lines: " << finalPolyData->GetNumberOfLines() << "\n"
	//	<< "  Verts: " << finalPolyData->GetNumberOfVerts() << "\n"
	//	<< "  Strips: " << finalPolyData->GetNumberOfStrips() << std::endl;
};

//void ScaffoldGenerator::print_cell_types() {
//	if (globalFaces[idx].size() != 3) {
//		std::cerr << "Non-triangle face found! Face #" << idx
//			<< " has " << globalFaces[idx].size() << " vertices." << std::endl;
//			<< " has " << globalFaces[idx].size() << " vertices." << std::endl;
//	}
//	else {
//		cells->InsertNextCell(3);
//		cells->InsertCellPoint(globalFaces[idx][0]);
//		cells->InsertCellPoint(globalFaces[idx][1]);
//		cells->InsertCellPoint(globalFaces[idx][2]);
//	}
//};

void ScaffoldGenerator::process_faces() {

	voro::c_loop_all cla(*con);
	voro::voronoicell_neighbor cell;

	int cellNr{ 1 };

	auto start_time = std::chrono::steady_clock::now();

	// Use a set of pairs to track unique faces by the IDs they separate
	std::cout << "Phase 1: Collecting unique face jobs..." << std::endl;
	std::set<std::pair<int, int>> processed_face_keys;

	if (cla.start()) do if (con->compute_cell(cell, cla)) {

		cellNr += 1;

		// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
		std::vector<double> cellVertices;
		std::vector<int> faceIndices;
		std::vector<int> cellNeighs;

		int seedId = cla.pid();

		// get position of seed and store it to an array
		double px = 0.0, py = 0.0, pz = 0.0;

		cla.pos(px, py, pz);

		// cell vertices in global system
		cell.vertices(px, py, pz, cellVertices);

		// nr of cell vertices
		int cellVertexNr = cellVertices.size() / 3;

		// get cell faces and neighbors
		cell.face_vertices(faceIndices);
		cell.neighbors(cellNeighs);

		// start processing faces 
		// loop in faces 
		int idx{ 0 };
		for (int i = 0; i < cellNeighs.size(); i++) {

			int order = faceIndices[idx];
			int neighborId = cellNeighs[i];

			if (neighborId < 0) { // Skip boundary faces
				idx += order + 1;
				continue;
			}

			// Create a canonical key from the two seed IDs
			std::pair<int, int> key = { std::min(seedId, neighborId), std::max(seedId, neighborId) };

			if (processed_face_keys.find(key) == processed_face_keys.end()) {
				processed_face_keys.insert(key);

				CellFace faceJob;

				// Compare factor (w) to the global user threshold 
				auto edgeData = graph->get_edge_width(seedId, neighborId);
				//float faceWidth = graph->get_edge_width(seedId, neighborId);
				
				// for edges with weight -1.0f
				if (edgeData.dstar < 0.0f) {
					std::cerr << "Warning: Edge not found in graph for face between "
						<< seedId << " and " << neighborId << ". Assuming solid wall." << std::endl;
					faceJob.createHole = false;
				}
				else {
					constexpr float eps = 1e-6f;
					faceJob.createHole = (edgeData.dstar >= thickness * 0.5 - 1e-6f);
					faceJob.delta = edgeData.delta;
				}

				for (int j = idx + 1; j < order + idx + 1; j++) {
					int localVertIdx = faceIndices[j];
					faceJob.localFace.push_back(localVertIdx);
					faceJob.vertices.push_back(cellVertices[3 * localVertIdx]);
					faceJob.vertices.push_back(cellVertices[3 * localVertIdx + 1]);
					faceJob.vertices.push_back(cellVertices[3 * localVertIdx + 2]);
				}
				allCellFaces.push_back(faceJob);
			}
			idx += order + 1;
		}
	} while (cla.inc());

	std::cout << "Phase 1 complete. Found " << allCellFaces.size() << " unique faces." << std::endl;

	// ======================================================================
	// PHASE 2: PROCESS FACES (Parallel)
	// ======================================================================
	std::cout << "Phase 2: Processing faces in parallel..." << std::endl;
	// Prepare output container for results (we'll collect results thread-locally and then merge)
	std::vector<NewFace> newFaces;
	newFaces.resize(allCellFaces.size()); // we optionally store per-index results

	// Use an auxiliary flag vector to mark successful faces (so we can skip empties later)
	std::vector<char> faceValid(allCellFaces.size(), 0);

	// Optional: tune number of threads via OMP_NUM_THREADS or omp_set_num_threads()
	#pragma omp parallel for 
	//private(triIdxs, holeIdxs, nodes, verts2D, tempVerts, holeVertices, chaikinVerts, finalVerts, finalVerts3d, u, v, origin, holeDll)
	for (int idxFace = 0; idxFace < static_cast<int>(allCellFaces.size()); ++idxFace) {

		const auto& cellFace = allCellFaces[idxFace];
		bool createHole = cellFace.createHole;

		// These variables are now automatically private to each thread
		std::vector<int> triIdxs;
		std::vector<int> holeIdxs;
		std::vector<Node*> nodes;
		Eigen::MatrixXd verts2D, tempVerts, holeVertices, chaikinVerts, finalVerts, finalVerts3d;
		Eigen::Vector3d u, v, origin;
		Cdll holeDll;

		// Minimal sanity: need at least 3 vertices (3 coords per vertex)
		if (cellFace.vertices.size() < 9) {
			// degenerate face
			continue;
		}

		std::vector<int> faceIndices = cellFace.localFace;
		Eigen::MatrixXd verts3D = vector_to_matrix3D(cellFace.vertices);
		ensure_ccw(faceIndices, verts3D);

		project_vertices_on_plane(verts3D, origin, u, v, verts2D);
		//project_vertices_on_plane(cellFace.vertices, origin, u, v, verts2D);

		if (verts2D.rows() < 3) continue; // still degenerate after projection

		// get the half spaces
		auto hs = get_poly_half_spaces(verts2D);

		// lambda function to estimate the left (inward) normal vectors
		auto rot90 = [](const Eigen::Vector2d& v) { return Eigen::Vector2d(-v.y(), v.x()); };

		// now based on threshold decide if we will drill the hole or not
		std::vector<std::vector<int>> faceCells;

		if (createHole){

			const double sfGlobal = pullbackRatio * (thickness / 2.0);

			holeVertices.resize(verts2D.rows(), 2);
			
			Eigen::Vector2d centroid = verts2D.colwise().mean();

			for (int vIdx = 0; vIdx < verts2D.rows(); ++vIdx) {
				int idx1 = (vIdx - 1 + verts2D.rows()) % verts2D.rows();
				int idx2 = vIdx;
				int idx3 = (vIdx + 1) % verts2D.rows();

				Eigen::Vector2d n1 = verts2D.row(idx2) - verts2D.row(idx1);
				Eigen::Vector2d n2 = verts2D.row(idx3) - verts2D.row(idx2);
				
				if (n1.norm() == 0 || n2.norm() == 0) {
					holeVertices.row(vIdx) = verts2D.row(vIdx);
					continue;
				}

				n1.normalize();
				n2.normalize();

				Eigen::Vector2d t1 = rot90(n1);
				Eigen::Vector2d t2 = rot90(n2);

				// current vertex
				Eigen::Vector2d xi = verts2D.row(idx2);

				// ensure they are inward the origin should ly above the half space
				if (t1.dot(centroid) > t1.dot(xi)) t1 = -t1;
				if (t2.dot(centroid) > t2.dot(xi)) t2 = -t2;

				//Eigen::Vector2d normal = t1 + t2;
				//if (normal.norm() > 1e-12) normal = t2;
				//else normal.normalize();

				Eigen::Vector2d normal = (t1.normalized() + t2.normalized()).normalized();

				//const double capLocal = 0.5 * std::min(n1.norm(), n2.norm());
				//const double sfLocal = std::min(sfGlobal, capLocal);

				//holeVertices(vIdx, 0) = verts2D(vIdx, 0) - sfLocal * normal.x();
				//holeVertices(vIdx, 1) = verts2D(vIdx, 1) - sfLocal * normal.y();

				double smax = std::numeric_limits<double>::infinity();
				for (const auto& L : hs) {
					double num = L.d - L.n.dot(xi);  // distance to line along its normal (>=0 inside)
					double den = L.n.dot(normal);    // how much moving along b approaches that line
					if (den > 1e-12) smax = std::min(smax, num / den);
				}
				if (!std::isfinite(smax)) smax = 0.0;

				double delta = cellFace.delta;
				double beta = pullbackRatio;
				double jitter = 0.0;
				double eta = (jitter > 0.0) ? (jitter * ((double)rand() / RAND_MAX * 2.0 - 1.0)) : 0.0;
				double si = std::clamp(delta * beta * (1.0 + eta), 0.0, smax);

				holeVertices.row(vIdx) = xi - 0.1 * normal;
			}

			if (holeVertices.rows() < 3 || holeVertices.rowwise().norm().maxCoeff() < 1e-6) {
				std::cout << " skipping " << std::endl;
				continue;
			}

			//Eigen::MatrixXd chaikinVerts;
			chaikin_subdivision(holeVertices, chaikinVerts);
			if (chaikinVerts.rows() < 3) {
				// fallback: use holeVertices directly
				chaikinVerts = holeVertices;
				if (chaikinVerts.rows() < 3) continue;
			}

			// find connector pair (closest point between outer polygon and inner polygon)
			Eigen::Vector2d pair = { 0, 0 };
			double minDist = std::numeric_limits<double>::infinity();
			for (int polyIdx = 0; polyIdx < verts2D.rows(); ++polyIdx) {
				for (int holeIdx = 0; holeIdx < chaikinVerts.rows(); ++holeIdx) {
					double dist = (chaikinVerts.row(holeIdx) - verts2D.row(polyIdx)).norm();
					if (dist < minDist) {
						pair = { static_cast<double>(polyIdx), static_cast<double>(holeIdx + verts2D.rows()) };
						minDist = dist;
					}
				}
			}
			if (minDist == std::numeric_limits<double>::infinity()) continue; // safety

			// build finalVerts (outer + hole points)
			finalVerts.resize(verts2D.rows() + chaikinVerts.rows(), 2);
			for (int cIdx = 0; cIdx < verts2D.rows(); ++cIdx) finalVerts.row(cIdx) = verts2D.row(cIdx);
			for (int cIdx = 0; cIdx < chaikinVerts.rows(); ++cIdx) finalVerts.row(verts2D.rows() + cIdx) = chaikinVerts.row(cIdx);

			// build a DLL of the hole indices and rotate to match connector - you already have this code
			nodes.clear();
			int startIdx = static_cast<int>(chaikinVerts.rows() + verts2D.rows()) - 1;
			for (int vIdx = startIdx; vIdx > verts2D.rows() - 1; --vIdx) {
				Node* ni = holeDll.append(vIdx);
				nodes.push_back(ni);
			}
			holeDll.rotate(chaikinVerts.rows() + verts2D.rows() - static_cast<int>(pair[1]) - 1);
			holeIdxs.clear();
			holeDll.get_data(holeIdxs);

			// Build triIdxs (indices to pass to ear clipping)
			triIdxs.clear();
			for (int cIdx = 0; cIdx < static_cast<int>(pair[0]) + 1; ++cIdx) triIdxs.push_back(cIdx);
			for (auto hIdx : holeIdxs) triIdxs.push_back(hIdx);
			for (int cIdx = static_cast<int>(pair[0]); cIdx < verts2D.rows(); ++cIdx) triIdxs.push_back(cIdx);

			// Ear clipping
			bool succed = ear_clipping(finalVerts, triIdxs, faceCells);
			if (!succed || faceCells.empty()) {
				// skip bad face
				continue;
			}

			// bring back to 3D
			finalVerts3d.resize(finalVerts.rows(), 3);
			back_to_3d(finalVerts3d, finalVerts, origin, u, v);

			// store into result slot (one per face)
			NewFace nf;
			nf.faces = std::move(faceCells);
			nf.vertices = std::move(finalVerts3d);

			newFaces[idxFace] = std::move(nf);
			faceValid[idxFace] = 1;

		}

		// if not selected to create a hole just do triangulation
		else {
			triIdxs.clear();
			for (int i = 0; i < verts2D.rows(); ++i) {
				triIdxs.push_back(i);
			}

			// Ear clipping on the simple (solid) polygon
			bool succed = ear_clipping(verts2D, triIdxs, faceCells);
			if (!succed || faceCells.empty()) {
				continue; // skip bad face
			}

			// bring back to 3D
			finalVerts3d.resize(verts2D.rows(), 3);
			back_to_3d(finalVerts3d, verts2D, origin, u, v);

			// store into result slot
			NewFace nf;
			nf.faces = std::move(faceCells);
			nf.vertices = std::move(finalVerts3d); // Note: using finalVerts3d
			newFaces[idxFace] = std::move(nf);
			faceValid[idxFace] = 1;
		}		
	}

	// After parallel region, compact results into a dense vector of NewFace for assembly
	std::vector<NewFace> compactFaces;
	compactFaces.reserve(allCellFaces.size());
	for (size_t i = 0; i < newFaces.size(); ++i) {
		if (faceValid[i]) compactFaces.push_back(std::move(newFaces[i]));
	}
	newFaces.swap(compactFaces);

	// --------------------------
	// Phase 3: assemble final globalVertices and globalFaces without locks
	// --------------------------
	std::cout << "Phase 3: Assembling final mesh..." << std::endl;

	// This check is important to remove any faces that failed during triangulation
	std::vector<NewFace> validFaces;
	for (auto& face : newFaces) {
		if (face.vertices.rows() > 0 && !face.faces.empty()) {
			validFaces.push_back(std::move(face));
		}
	}

	// Clear old data and build the final mesh from scratch
	globalVertices.clear();
	globalFaces.clear();
	globalVertexMap.clear();
	globalIndex = 0;

	// Loop through the valid results and stitch them together
	for (const auto& face : validFaces) {

		// This map will store the mapping from this face's local vertex indices
		// to their new global indices.
		std::vector<int> current_face_global_indices;
		current_face_global_indices.resize(face.vertices.rows());

		// First, process the vertices for this face, checking for duplicates
		for (int i = 0; i < face.vertices.rows(); ++i) {
			std::array<double, 3> coords = { face.vertices(i, 0), face.vertices(i, 1), face.vertices(i, 2) };
			GlobalVertex key(coords);

			auto it = globalVertexMap.find(key);
			if (it == globalVertexMap.end()) {
				// This is a new vertex
				globalVertexMap[key] = globalIndex;
				current_face_global_indices[i] = globalIndex;
				globalVertices.push_back(coords);
				globalIndex++;
			}
			else {
				// This vertex already exists, reuse its index
				current_face_global_indices[i] = it->second;
			}
		}

		// Now, create the global faces using the re-mapped indices
		for (const auto& local_triangle : face.faces) {
			globalFaces.push_back({
				current_face_global_indices[local_triangle[0]],
				current_face_global_indices[local_triangle[1]],
				current_face_global_indices[local_triangle[2]]
				});
		}
	}

	std::cout << "Assembled mesh: vertices=" << globalVertices.size()
		<< ", triangles=" << globalFaces.size() << std::endl;

	delete con;
	//render_vtk_face(globalVertices, globalFaces, "FinalMesh");

	Eigen::MatrixXd meshVertices(globalVertices.size(), 3);
	for (int vrow{ 0 }; vrow < globalVertices.size(); vrow++) {
		Eigen::Vector3d row = { globalVertices[vrow][0], globalVertices[vrow][1], globalVertices[vrow][2] };
		meshVertices.row(vrow) = row;
	}

	// Create points
	vtkNew<vtkPoints> points;

	if (meshVertices.cols() == 3) {
		for (int i = 0; i < meshVertices.rows(); ++i) {
			points->InsertNextPoint(meshVertices(i, 0), meshVertices(i, 1), meshVertices(i, 2));
		}
	}
	else if (meshVertices.cols() == 2) {
		for (int i = 0; i < meshVertices.rows(); ++i) {
			points->InsertNextPoint(meshVertices(i, 0), meshVertices(i, 1), 0.0);
		}
	}

	vtkNew<vtkCellArray> cells;

	for (int idx{ 0 }; idx < globalFaces.size(); idx++) {
		if (globalFaces[idx].size() != 3) {
			std::cerr << "Non-triangle face found! Face #" << idx
				<< " has " << globalFaces[idx].size() << " vertices." << std::endl;
		}
		else {
			cells->InsertNextCell(3);
			cells->InsertCellPoint(globalFaces[idx][0]);
			cells->InsertCellPoint(globalFaces[idx][1]);
			cells->InsertCellPoint(globalFaces[idx][2]);
		}
	}

	// Update scaffold polydata object
	scaffoldMesh = vtkSmartPointer<vtkPolyData>::New();
	scaffoldMesh->SetPoints(points);
	scaffoldMesh->SetPolys(cells);

	scaffoldMesh->BuildCells();
	scaffoldMesh->BuildLinks();

	auto end_time = std::chrono::steady_clock::now();

	// Calculate the duration in milliseconds
	auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

	// Format the output string
	std::cout << " face processing ended in "
		<< std::fixed << std::setprecision(3) // Set precision to 3 decimal places
		<< duration_ms.count() / 1000.0   // Convert ms to seconds
		<< " seconds." << std::endl;

};

std::vector<unsigned int> ScaffoldGenerator::get_connected_edges() {
	return connectedIndices;
};

std::vector<float> ScaffoldGenerator::get_connected_vertices() {
	return connectedVertices;
};

// -------------------------------------------------------------------------------------
// Scaffold Generator inside a rectangle
ScaffoldGeneratorBox::ScaffoldGeneratorBox(
	std::vector<std::array<double, 3>>& seeds,
	const std::array<float, 6>& bounds,
	const std::array<int, 3>& blockDim,
	const double thickness,
	const double pullbackRatio,
	std::function<bool(const std::array<double, 3>&)> is_inside
) : ScaffoldGenerator(seeds, bounds, blockDim, thickness, pullbackRatio, is_inside) {

	con = new voro::container(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
		blockDim[0], blockDim[1], blockDim[2], false, false, false, 16
	);

};

void ScaffoldGenerator::regularize_voro(int regSteps) {

	std::cout << "Regularization Steps: " << regSteps << std::endl;

	for (int step{ 0 }; step < regSteps; ++step) {

		std::cout << "Regularization Step: " << step + 1 << std::endl;

		std::vector<std::array<double, 3>> newSeeds;

		// clear previous container
		con->clear();

		for (int i{ 0 }; i < seeds.size(); i++) {
			con->put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
		}

		voro::c_loop_all cla(*con);
		voro::voronoicell_neighbor cell;

		if (cla.start()) do if (con->compute_cell(cell, cla)) {

			//std::cout << "Entered" << std::endl;
			// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
			std::vector<double> cellVertices;
			std::vector<int> faceVertices;
			std::vector<int> cellNeighs;

			int seedId = cla.pid();

			// get position of seed and store it to an array
			double px = 0.0, py = 0.0, pz = 0.0;

			cla.pos(px, py, pz);

			// --------------------

			// cell vertices in global system
			cell.vertices(px, py, pz, cellVertices);

			// get cell faces and neighbors
			cell.face_vertices(faceVertices);
			cell.neighbors(cellNeighs);

			std::array<double, 3> centroid;
			cell.centroid(centroid[0], centroid[1], centroid[2]);

			centroid[0] += px;
			centroid[1] += py;
			centroid[2] += pz;
			
			std::cout << "here" << std::endl;
			if (is_inside({ centroid[0], centroid[1], centroid[2] })) {
				newSeeds.push_back(centroid);
			}

		} while (cla.inc());

		seeds = newSeeds;

	}

	std::cout << "------------------" << std::endl;
	std::cout << "after reg seeds: " << seeds.size() << std::endl;
}

void ScaffoldGeneratorBox::populate_voro(const int regSteps) {


	con->clear();

	std::vector<vtkSmartPointer<vtkPolyData>> polys;

	for (int i{ 0 }; i < seeds.size(); i++) {
		con->put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
	}

	if (regSteps > 1) {
		regularize_voro(regSteps);
	}

	// populate the graph
	voro::c_loop_all cla(*con);
	voro::voronoicell_neighbor cell;

	graph = std::make_unique<Graph>();

	// populat the graph nodes
	if (cla.start()) do {
		graph->add_vertex(cla.pid());
	} while (cla.inc());

	// loop through the cells
	if (cla.start()) do if (con->compute_cell(cell, cla)) {

		// get the id of the seed
		int seedId = cla.pid();

		//double x{ 0.0 }, y{ 0.0 }, z{ 0.0 };

		//cla.pos(x, y, z);

		// get the neighbors and the corresponding face area
		std::vector<int> cellNeighs;
		cell.neighbors(cellNeighs);

		std::vector<int> faceIndices;
		cell.face_vertices(faceIndices);

		double px{ 0.0 }, py{ 0.0 }, pz{ 0.0 };
		std::vector<double> cellVertices;
		cla.pos(px, py, pz);
		cell.vertices(px, py, pz, cellVertices);

		// get each face areas
		std::vector<double> faceAreas;
		cell.face_areas(faceAreas);

		// get each face perimeter
		//std::vector<double> facePerimeters;
		//cell.face_perimeters(facePerimeters);
		std::vector<double> faceNormalsRaw;
		cell.normals(faceNormalsRaw);

		const int numFaces = static_cast<int>(cellNeighs.size());
		std::vector<double> faceMinWidths(numFaces, 0.0);
		//std::vector<double> faceMaxWidths(numFaces, 0.0);
		std::vector<Eigen::Vector3d> faceNormals(numFaces, Eigen::Vector3d::Zero());
		std::vector<double> planeCoeffs(numFaces, 0.0); 
		std::vector<Eigen::Vector3d> faceCentroids(numFaces, Eigen::Vector3d::Zero());

		int idxOffset = 0;

		// --- Pass 1: Compute minWidth and normals for each face ---
		for (int i = 0; i < numFaces; ++i) {

			int neighborId = cellNeighs[i];
			int order = faceIndices[idxOffset];

			// avoid process of cells on the boundaries
			//if (neighborId >= 0) {

				std::vector<double> faceVertices3D;		
				
				for (int j = 1; j <= order; ++j) {
					int localVertIdx = faceIndices[idxOffset + j];
					faceVertices3D.push_back(cellVertices[3 * localVertIdx + 0]);
					faceVertices3D.push_back(cellVertices[3 * localVertIdx + 1]);
					faceVertices3D.push_back(cellVertices[3 * localVertIdx + 2]);
				}

				if (3 * i + 2 < static_cast<int>(faceNormalsRaw.size())) {
					faceNormals[i] = Eigen::Vector3d(
						faceNormalsRaw[3 * i + 0],
						faceNormalsRaw[3 * i + 1],
						faceNormalsRaw[3 * i + 2]
					).normalized();
				}
				else {
					faceNormals[i] = Eigen::Vector3d::Zero();
				}

				// add the to connected vertices also the mean (centroid) of the face only once
				Eigen::Vector3d centroid = vector_to_matrix3D(faceVertices3D).colwise().mean().transpose();
				faceCentroids[i] = centroid;

				if (neighborId >= 0 && seedId < neighborId) {
					centroids[seedId].insert({ neighborId, centroid });
				}

				// estimate the plane coefficient d
				planeCoeffs[i] = faceNormals[i].dot(centroid);
	
				//Eigen::MatrixXd verts2D;
				//Eigen::Vector3d u, v, origin;
				//project_vertices_on_plane(faceVertices3D, origin, u, v, verts2D);

				//if (verts2D.rows() >= 3) {
				//	//double minWidth = polygon_min_width(verts2D);

				//	PolygonWidth w = get_polygon_width(verts2D);

				//	double minWidth = w.minWidth;
				//	double maxWidth = w.maxWidth;
				//	faceMinWidths[i] = minWidth;
				//	//faceMaxWidths[i] = minWidth;
				//	graph->add_edge(seedId, neighborId, static_cast<float>(minWidth));

				//}


			//}

			idxOffset += order + 1;
		}
		
		// --- Pass 2: For each face check if the eroded cell with block it or not
		idxOffset = 0;
		for (int i = 0; i < numFaces; ++i) {

			int neighborId = cellNeighs[i];
			int order = faceIndices[idxOffset];
			
			//if (neighborId >= 0) {
				
				Eigen::Vector3d normal = faceNormals[i];

				std::vector<double> faceVertices3D;

				for (int j = 1; j <= order; ++j) {
					int localVertIdx = faceIndices[idxOffset + j];
					faceVertices3D.push_back(cellVertices[3 * localVertIdx + 0]);
					faceVertices3D.push_back(cellVertices[3 * localVertIdx + 1]);
					faceVertices3D.push_back(cellVertices[3 * localVertIdx + 2]);
				}

				Eigen::MatrixXd verts2D; Eigen::Vector3d u, v, origin;
				project_vertices_on_plane(faceVertices3D, origin, u, v, verts2D);

				//double dstar = survival_margin(faceNormals, planeCoeffs, ring3D, i);

				////double dstar = polygon_is_open(
				//bool isOpen = polygon_is_open(
				//	faceNormals,
				//	planeCoeffs,
				//	verts2D,
				//	i,
				//	thickness * 0.5
				//);


				//std::cout << isOpen << std::endl;
				//double weight = (isOpen) ? 1.0 : 0.0;

				Eigen::Matrix<double, 3, 2> basis;
				basis.col(0) = u;
				basis.col(1) = v;

				// estimate the delta as the inlet radius
				double delta =  pullbackRatio * polygon_inradius(verts2D);

				double d = erosion_margin_for_face(
					i, 
					faceNormals,
					planeCoeffs,
					basis,
					origin,
					verts2D,
					delta
				);

				//std::cout << d << " " << delta << std::endl;
				
				//bool open = dstar > 0.0 + 1e-6;
				graph->add_edge(seedId, neighborId,
					{ static_cast<float>(d), static_cast<float>(delta) });
			//}
			
			idxOffset += order + 1;
		}

	} while (cla.inc());
};

void ScaffoldGeneratorBox::add_cylindrical_wall(
	const double pt0, const double pt1, const double pt2,
	const double axis0,
	const double axis1,
	const double axis2,
	const double radius
) {
	wall = std::make_unique<voro::wall_cylinder>(pt0, pt1, pt2, axis0, axis1, axis2, radius);
	con->add_wall(*wall);
	std::cout << "added wall" << std::endl;
}

void ScaffoldGeneratorBox::check_global_vertices(const Eigen::MatrixXd& vertices, std::vector<int>* face) {

	for (int rowIdx{ 0 }; rowIdx < vertices.rows(); rowIdx++) {
		std::array<double, 3> coords = {
			vertices(rowIdx, 0),
			vertices(rowIdx, 1),
			vertices(rowIdx, 2),
		};

		GlobalVertex key(coords);
		auto it = globalVertexMap.find(key);
		if (it == globalVertexMap.end()) {
			globalVertexMap[key] = globalIndex;
			localToGlobal[rowIdx] = globalIndex;
			if (face) {
				face->push_back(globalIndex);
			}
			globalVertices.push_back({ coords[0], coords[1], coords[2] });
			globalIndex++;
		}
		else {
			localToGlobal[rowIdx] = it->second;
			if (face) {
				face->push_back(it->second);
			}
		}
	}
};

// -----------------------------------------------------------------------------------------------------
// generate mesh inside wall

ScaffoldGeneratorWall::ScaffoldGeneratorWall(
	std::vector<std::array<double, 3>>& seeds,
	vtkSmartPointer<vtkPolyData>& containerPoly,
	const std::array<int, 3>& blockDim, 
	const int neighbors,
	const float minDist,
	const double thickness,
	const double pullbackRatio,
	std::function<bool(const std::array<double, 3>&)> is_inside
	) : ScaffoldGenerator(seeds, { 0,0,0,0,0,0 }, blockDim, thickness, pullbackRatio, is_inside), containerMesh(containerPoly), neighbors(neighbors), minDist(minDist) {

	double bds[6];

	double cap = 0.5;

	containerMesh->GetBounds(bds);
	bounds[0] = bds[0] - cap;
	bounds[1] = bds[1] + cap;
	bounds[2] = bds[2] - cap;
	bounds[3] = bds[3] + cap;
	bounds[4] = bds[4] - cap;
	bounds[5] = bds[5] + cap;

	_process_triangles();
};

void ScaffoldGeneratorWall::generate_voro(const int regSteps) {

	bool kdsurface{ false };

	std::cout << "Seed Nr: " << seeds.size() << std::endl;
	std::cout << "Barycenter Nr: " << bCenters.size() << std::endl;
	std::cout << "Normals Nr: " << normals.size() << std::endl;

	// create kdtree to locate points 
	// Create and build the k-d tree
	vtkSmartPointer<vtkKdTreePointLocator> kdTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	if (kdsurface) {
		kdTree->SetDataSet(containerMesh);
	}
	else {
		vtkSmartPointer<vtkPoints> bcvtk = vtkSmartPointer<vtkPoints>::New();
		for (int i{ 0 }; i < bCenters.size(); i++) {
			bcvtk->InsertNextPoint(bCenters[i][0], bCenters[i][1], bCenters[i][2]);
		}

		vtkNew<vtkPolyData> bcvtkdata;
		bcvtkdata->SetPoints(bcvtk);
		kdTree->SetDataSet(bcvtkdata);
	}
	kdTree->BuildLocator();

	// create container
	con = new voro::container(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
		blockDim[0], blockDim[1], blockDim[2], false, false, false, 16);

	MeshWall wall(containerMesh, normals, kdTree, neighbors);
	con->add_wall(wall);

	std::vector<vtkSmartPointer<vtkPolyData>> polys;

	//regularize_voro(regSteps, containerMesh);

	//// clear previous container
	//con->clear();

	for (int i{ 0 }; i < seeds.size(); i++) {
		con->put(i, seeds[i][0], seeds[i][1], seeds[i][2]);
	}

	//for (int i{ 0 }; i < bCenters.size(); i++) {
	//	seeds.push_back({ bCenters[i][0], bCenters[i][1], bCenters[i][2] });
	//	con->put(i, bCenters[i][0], bCenters[i][1], bCenters[i][2]);
	//}
	int originalSeedCount = seeds.size();  // store original seed count

	for (int i = 0; i < bCenters.size(); i++) {
		seeds.push_back({ bCenters[i][0], bCenters[i][1], bCenters[i][2] });
		con->put(originalSeedCount + i, bCenters[i][0], bCenters[i][1], bCenters[i][2]);
	}

	//std::cout << "Nr of seeds: " << seeds.size() << std::endl;
	//process_faces(scaleFactor, edgeSize);
}

void ScaffoldGeneratorWall::check_global_vertices(const Eigen::MatrixXd& vertices, std::vector<int>* face) {

	for (int rowIdx{ 0 }; rowIdx < vertices.rows(); rowIdx++) {
		std::array<double, 3> coords = {
			vertices(rowIdx, 0),
			vertices(rowIdx, 1),
			vertices(rowIdx, 2),
		};
		GlobalVertex key(coords);
		auto it = globalVertexMap.find(key);
		if (it == globalVertexMap.end()) {
			globalVertexMap[key] = globalIndex;
			localToGlobal[rowIdx] = globalIndex;
			if (face) {
				face->push_back(globalIndex);
			}
			globalVertices.push_back({ coords[0], coords[1], coords[2] });
			globalIndex++;
		}
		else {
			localToGlobal[rowIdx] = it->second;
			if (face) {
				face->push_back(it->second);
			}
		}
	}
};


void ScaffoldGeneratorWall::_process_triangles() {

	std::cout << "processing triangles" << std::endl;
	
	std::vector<std::array<double, 3>> centers;
	std::vector<std::array<double, 3>> tempNormals;
	for (vtkIdType i = 0; i < containerMesh->GetNumberOfCells(); ++i) {

		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(containerMesh->GetCell(i));
		if (!triangle) continue;

		double p1[3], p2[3], p3[3];
		triangle->GetPoints()->GetPoint(0, p1);
		triangle->GetPoints()->GetPoint(1, p2);
		triangle->GetPoints()->GetPoint(2, p3);

		//trianglePt.push_back({p1[0], p1[1], p1[2]});

		// convert pts to Eigen
		Eigen::Vector3d pe1{ p1[0], p1[1], p1[2] };
		Eigen::Vector3d pe2{ p2[0], p2[1], p2[2] };
		Eigen::Vector3d pe3{ p3[0], p3[1], p3[2] };

		// estimate the barycenter
		double bc1 = (p1[0] + p2[0] + p3[0]) / 3;
		double bc2 = (p1[1] + p2[1] + p3[1]) / 3;
		double bc3 = (p1[2] + p2[2] + p3[2]) / 3;
		std::array<double, 3> bc = { bc1, bc2, bc3 };

		Eigen::Vector3d normal = (pe2 - pe1).cross(pe3 - pe1);
		normal.normalize();
		tempNormals.push_back({ normal[0], normal[1], normal[2] });

		bc[0] += 0.5 * normal[0];
		bc[1] += 0.5 * normal[1];
		bc[2] += 0.5 * normal[2];

		centers.push_back(bc);
		//seeds.push_back(bc);
	}

	// filter centers based on the minimum distance

	// Filtering step
	for (size_t i = 0; i < centers.size(); i++) {
		bool keep = true;
		for (size_t j = 0; j < bCenters.size(); j++) {
			if (distance(centers[i], bCenters[j]) < minDist) {
				keep = false;
				break;
			}
		}
		if (keep) {
			bCenters.push_back(centers[i]);
			normals.push_back(tempNormals[i]);
		}
	}
	std::cout << "Before filtering: " << centers.size() << std::endl;
	std::cout << "After filtering: " << bCenters.size() << std::endl;
};


// ------------------------------------------------------------------
// Volume Optimization

VolOpt::VolOpt(
	std::vector<std::array<double, 3>>& seeds,
	const Eigen::VectorXd targets,
	const Eigen::VectorXd initV,
	const std::array<double, 6>& bounds,
	std::function<void(const std::string&)> logCallback) :
		currSeeds(seeds), targetVols(targets), wInit(initV), bounds(bounds), log_callback(logCallback) {};

void VolOpt::loop(const int regSteps) {

	for (int step = 0; step < regSteps; step++) {

		if (log_callback) {
			log_callback("Regularization step: " + std::to_string(step + 1));
		}
		else {
			std::cout << "Regularization step : " << step + 1 << std::endl;
		}

		// create an instance of the objective function
		myFunc func(
			bounds,
			currSeeds,
			wInit,
			targetVols
		);

		double f0;
		Eigen::VectorXd g0 = Eigen::VectorXd::Zero(wInit.size());
		f0 = func(wInit, g0);

		// estimate the tolerance based on the initial step gradient
		double tol = std::min(0.1 * 1e-2 * targetVols.minCoeff(), 0.1 * 1e-2 * targetVols.minCoeff() / g0.array().abs().maxCoeff());

		if (log_callback) {
			log_callback("Optimality Tolerance: " + std::to_string(tol));
		}
		else {
			std::cout << "Optimality Tolerance: " << tol << std::endl;
		}

		// create a solver object
		BFGS<myFunc> solver(func);

		solver.maxIter = 1000;
		//solver.epsilon = 0.00056;
		solver.epsilon = tol;
		solver.minimize(wInit);

		//std::cout << func.w << std::endl;
		//log_callback("" + std::to_string(func.volError));
		std::cout << "Volume Error: " << func.volError << std::endl;
		//std::cout << func.w.transpose() << std::endl;

		// get the final cells
		// Convert weights to radii
		Eigen::VectorXd radii = convert_radii(func.w);

		std::cout << "Creating Container" << std::endl;
		// create container
		voro::container_poly* conp = new voro::container_poly(
			bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
			10, 10, 10,
			false, false, false, 16
		);

		for (int i = 0; i < currSeeds.size(); i++) {
			conp->put(i, func.centroids[i][0], func.centroids[i][1], func.centroids[i][2], radii[i]);
		}

		Eigen::VectorXd currVols = Eigen::VectorXd::Zero(currSeeds.size());
		Eigen::VectorXd currCosts = Eigen::VectorXd::Zero(currSeeds.size());

		// loop in cells
		voro::c_loop_all cla(*conp);
		voro::voronoicell_neighbor cell;

		if (cla.start()) do if (conp->compute_cell(cell, cla)) {

			// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
			std::vector<double> cellVertices;
			std::vector<int> faceVertices;
			std::vector<int> cellNeighs;

			int seedId = cla.pid();

			// get position of seed and store it to an array
			double x = 0.0, y = 0.0, z = 0.0;
			cla.pos(x, y, z);

			// cell vertices in global system
			cell.vertices(x, y, z, cellVertices);

			// get cell faces and neighbors
			cell.face_vertices(faceVertices);
			cell.neighbors(cellNeighs);

			vtkSmartPointer<vtkPolyData> poly = cell_2_vtk(cellNeighs, cellVertices, faceVertices);
			polys.push_back(poly);

		} while (cla.inc());

		delete conp;

		// check if the centroids and the currSeeds are close enough
		double maxDiff{ 0 };
		for (int i{ 0 }; i < wInit.size(); i++) {
			double diff = std::sqrt(
				std::pow(currSeeds[i][0] - func.centroids[i][0], 2) +
				std::pow(currSeeds[i][1] - func.centroids[i][1], 2) +
				std::pow(currSeeds[i][2] - func.centroids[i][2], 2));

			if (maxDiff < diff) { maxDiff = diff; }
		}
		//std::cout << "max diff between centroids and previous seeds: " << maxDiff << std::endl;

		//log_callback("max diff between centroids and previous seeds: " + std::to_string(maxDiff));

		if (maxDiff < 0.001) {
			//log_callback("Regularization ending!");
			break;
		}
		else if (step < regSteps - 1) {
			polys.clear();
		}
		else {
			break;
		}

		// update seeds

		currSeeds = func.centroids;

		//log_callback("\n ----------- \n");
	}
	vtkNew<vtkAppendPolyData> appendFilter;
	for (int i = 0; i < polys.size(); i++) {
		appendFilter->AddInputData(polys[i]);
	}
	appendFilter->Update();
	//scaffoldMesh = vtkSmartPointer<vtkPolyData>::New();
	scaffoldMesh = appendFilter->GetOutput();
};

void VolOpt::generate_mesh(
	const double& thickness,
	vtkSmartPointer<vtkPolyData>& finalPolyData,
	const std::array<int, 3>& res) {

	vtkNew<vtkNamedColors> colors;

	// decide mesh resolution
	int dim[3] = {};

	if (res.empty()) {
		if (thickness < 0.3) {
			dim[0] = 300;
			dim[1] = 300;
			dim[2] = 300;
		}
		else if (0.3 <= thickness && thickness < 0.5) {
			dim[0] = 200;
			dim[1] = 200;
			dim[2] = 200;
		}
		else {
			dim[0] = 100;
			dim[1] = 100;
			dim[2] = 100;
		}
	}
	else {
		dim[0] = res[0];
		dim[1] = res[1];
		dim[2] = res[2];
	}
	// build along line
	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->AddInputData(scaffoldMesh);
	implictModeller->SetMaximumDistance(thickness / 2 * 1.01);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	implictModeller->SetModelBounds(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

	// extract isosurface
	vtkNew<vtkContourFilter> isoFilter;
	isoFilter->SetInputConnection(implictModeller->GetOutputPort());
	isoFilter->SetValue(0, thickness / 2.0f);
	isoFilter->Update();

	// update normals
	vtkNew<vtkPolyDataNormals> norms;
	//vtkPolyDataNormals* norms = vtkPolyDataNormals::New();

	norms->SetInputConnection(isoFilter->GetOutputPort());
	norms->ComputePointNormalsOff();
	norms->ComputeCellNormalsOn();
	norms->ConsistencyOn();
	norms->AutoOrientNormalsOn();
	norms->Update();

	// use window sinc filter
	//vtkNew<vtkWindowedSincPolyDataFilter> sincSmoother;
	//sincSmoother->SetInputConnection(norms->GetOutputPort());
	//sincSmoother->SetNumberOfIterations(30);
	//sincSmoother->BoundarySmoothingOn();
	//sincSmoother->FeatureEdgeSmoothingOff();
	//sincSmoother->Update();
	
	finalPolyData = norms->GetOutput();
	//vtkNew<vtkSTLWriter> stlWriter;
	//stlWriter->SetFileName(fileName.c_str());
	////stlWriter->SetInputConnection(sincSmoother->GetOutputPort());
	//stlWriter->SetInputConnection(norms->GetOutputPort());
	//stlWriter->Write();
};

void VolOpt::get_seeds(std::vector<std::array<double, 3>>& outSeeds) {
	outSeeds = currSeeds;
};

// ----------------------------------------------------------------
VolOptWall::VolOptWall(
	std::vector<std::array<double, 3>>& seeds,
	const Eigen::VectorXd targets,
	const Eigen::VectorXd initV,
	vtkSmartPointer<vtkPolyData>& containerPoly,
	int neighbors,
	std::function<void(const std::string&)> logCallback) :
	currSeeds(seeds), targetVols(targets), wInit(initV), containerMesh(containerPoly),
	neighbors(neighbors), log_callback(logCallback) {

	double bds[6];
	containerMesh->GetBounds(bds);

	bounds[0] = bds[0];
	bounds[1] = bds[1];
	bounds[2] = bds[2];
	bounds[3] = bds[3];
	bounds[4] = bds[4];
	bounds[5] = bds[5];

	_process_triangles();
};

void VolOptWall::loop(const int regSteps) {

	// create container
	bool kdsurface{ false };
	vtkSmartPointer<vtkKdTreePointLocator> kdTree = vtkSmartPointer<vtkKdTreePointLocator>::New();
	if (kdsurface) {
		kdTree->SetDataSet(containerMesh);
	}
	else {
		vtkSmartPointer<vtkPoints> bcvtk = vtkSmartPointer<vtkPoints>::New();
		for (int i{ 0 }; i < bCenters.size(); i++) {
			bcvtk->InsertNextPoint(bCenters[i][0], bCenters[i][1], bCenters[i][2]);
		}

		vtkNew<vtkPolyData> bcvtkdata;
		bcvtkdata->SetPoints(bcvtk);
		kdTree->SetDataSet(bcvtkdata);
	}
	kdTree->BuildLocator();

	MeshWall abstractWall(containerMesh, normals, kdTree, neighbors);

	for (int step = 0; step < regSteps; step++) {

		if (log_callback) {
			log_callback("Regularization step: " + std::to_string(step + 1));
		}
		else {
			std::cout << "Regularization step : " << step + 1 << std::endl;
		}

		// create an instance of the objective function
		myFuncWall func(
			abstractWall,
			bounds,
			currSeeds,
			wInit,
			targetVols
		);

		double f0;
		Eigen::VectorXd g0 = Eigen::VectorXd::Zero(wInit.size());
		f0 = func(wInit, g0);

		// estimate the tolerance based on the initial step gradient
		double tol = std::min(0.1 * 1e-2 * targetVols.minCoeff(), 0.1 * 1e-2 * targetVols.minCoeff() / g0.array().abs().maxCoeff());

		if (log_callback) {
			log_callback("Optimality Tolerance: " + std::to_string(tol));
		}
		else {
			std::cout << "Optimality Tolerance: " << tol << std::endl;
		}

		// create a solver object
		BFGS<myFuncWall> solver(func);

		solver.maxIter = 1000;
		//solver.epsilon = 0.00056;
		solver.epsilon = tol;
		solver.minimize(wInit);

		//std::cout << func.w << std::endl;
		//log_callback("" + std::to_string(func.volError));
		std::cout << "Volume Error: " << func.volError << std::endl;
		//std::cout << func.w.transpose() << std::endl;

		// get the final cells
		// Convert weights to radii
		Eigen::VectorXd radii = convert_radii(func.w);

		std::cout << "Creating Container" << std::endl;

		voro::container_poly* conp = new voro::container_poly(
			bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5],
			10, 10, 10,
			false, false, false, 16
		);

		for (int i = 0; i < currSeeds.size(); i++) {
			conp->put(i, func.centroids[i][0], func.centroids[i][1], func.centroids[i][2], radii[i]);
		}

		conp->add_wall(abstractWall);

		Eigen::VectorXd currVols = Eigen::VectorXd::Zero(currSeeds.size());
		Eigen::VectorXd currCosts = Eigen::VectorXd::Zero(currSeeds.size());

		// loop in cells
		voro::c_loop_all cla(*conp);
		voro::voronoicell_neighbor cell;

		if (cla.start()) do if (conp->compute_cell(cell, cla)) {

			// vector to store the cell vertexes in global coordinates, cell neighbors and face vertices
			std::vector<double> cellVertices;
			std::vector<int> faceVertices;
			std::vector<int> cellNeighs;

			int seedId = cla.pid();

			// get position of seed and store it to an array
			double x = 0.0, y = 0.0, z = 0.0;
			cla.pos(x, y, z);

			// cell vertices in global system
			cell.vertices(x, y, z, cellVertices);

			// get cell faces and neighbors
			cell.face_vertices(faceVertices);
			cell.neighbors(cellNeighs);

			vtkSmartPointer<vtkPolyData> poly = cell_2_vtk(cellNeighs, cellVertices, faceVertices);
			polys.push_back(poly);

		} while (cla.inc());

		delete conp;

		// check if the centroids and the currSeeds are close enough
		double maxDiff{ 0 };
		for (int i{ 0 }; i < wInit.size(); i++) {
			double diff = std::sqrt(
				std::pow(currSeeds[i][0] - func.centroids[i][0], 2) +
				std::pow(currSeeds[i][1] - func.centroids[i][1], 2) +
				std::pow(currSeeds[i][2] - func.centroids[i][2], 2));

			if (maxDiff < diff) { maxDiff = diff; }
		}
		//std::cout << "max diff between centroids and previous seeds: " << maxDiff << std::endl;

		//log_callback("max diff between centroids and previous seeds: " + std::to_string(maxDiff));

		if (maxDiff < 0.001) {
			//log_callback("Regularization ending!");
			break;
		}
		else if (step < regSteps - 1) {
			polys.clear();
		}
		else {
			break;
		}

		// update seeds

		currSeeds = func.centroids;

		//log_callback("\n ----------- \n");
	}
	vtkNew<vtkAppendPolyData> appendFilter;
	for (int i = 0; i < polys.size(); i++) {
		appendFilter->AddInputData(polys[i]);
	}
	appendFilter->Update();
	//scaffoldMesh = vtkSmartPointer<vtkPolyData>::New();
	scaffoldMesh = appendFilter->GetOutput();
};

void VolOptWall::_process_triangles() {

	std::cout << "processing triangles" << std::endl;

	for (vtkIdType i = 0; i < containerMesh->GetNumberOfCells(); ++i) {

		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(containerMesh->GetCell(i));
		if (!triangle) continue;

		double p1[3], p2[3], p3[3];
		triangle->GetPoints()->GetPoint(0, p1);
		triangle->GetPoints()->GetPoint(1, p2);
		triangle->GetPoints()->GetPoint(2, p3);

		//trianglePt.push_back({p1[0], p1[1], p1[2]});

		// convert pts to Eigen
		Eigen::Vector3d pe1{ p1[0], p1[1], p1[2] };
		Eigen::Vector3d pe2{ p2[0], p2[1], p2[2] };
		Eigen::Vector3d pe3{ p3[0], p3[1], p3[2] };

		// estimate the barycenter
		double bc1 = (p1[0] + p2[0] + p3[0]) / 3;
		double bc2 = (p1[1] + p2[1] + p3[1]) / 3;
		double bc3 = (p1[2] + p2[2] + p3[2]) / 3;
		std::array<double, 3> bc = { bc1, bc2, bc3 };

		bCenters.push_back(bc);
		currSeeds.push_back(bc);

		Eigen::Vector3d normal = (pe2 - pe1).cross(pe3 - pe1);
		normal.normalize();
		normals.push_back({ normal[0], normal[1], normal[2] });
	}
};

void VolOptWall::generate_mesh(
	const double& thickness,
	const std::string& fileName, 
	const std::array<int, 3>& res) {

	vtkNew<vtkNamedColors> colors;

	// decide mesh resolution
	int dim[3] = {};

	if (res.empty()) {
		if (thickness < 0.3) {
			dim[0] = 300;
			dim[1] = 300;
			dim[2] = 300;
		}
		else if (0.3 <= thickness && thickness < 0.5) {
			dim[0] = 100;
			dim[1] = 100;
			dim[2] = 100;
		}
		else {
			dim[0] = 100;
			dim[1] = 100;
			dim[2] = 100;
		}
	}
	else {
		dim[0] = res[0];
		dim[1] = res[1];
		dim[2] = res[2];
	}

	// build along line
	vtkNew<vtkImplicitModeller> implictModeller;
	implictModeller->AddInputData(scaffoldMesh);
	implictModeller->SetMaximumDistance(thickness / 2 * 1.01);
	implictModeller->SetSampleDimensions(dim[0], dim[1], dim[2]);
	implictModeller->SetModelBounds(
		bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);

	// extract isosurface
	vtkNew<vtkContourFilter> isoFilter;
	isoFilter->SetInputConnection(implictModeller->GetOutputPort());
	isoFilter->SetValue(0, thickness / 2);
	isoFilter->Update();

	// update normals
	vtkNew<vtkPolyDataNormals> norms;
	//vtkPolyDataNormals* norms = vtkPolyDataNormals::New();

	norms->SetInputConnection(isoFilter->GetOutputPort());
	norms->ComputePointNormalsOff();
	norms->ComputeCellNormalsOn();
	norms->ConsistencyOn();
	norms->AutoOrientNormalsOn();
	norms->Update();

	// use window sinc filter
	//vtkNew<vtkWindowedSincPolyDataFilter> sincSmoother;
	//sincSmoother->SetInputConnection(norms->GetOutputPort());
	//sincSmoother->SetNumberOfIterations(30);
	//sincSmoother->BoundarySmoothingOn();
	//sincSmoother->FeatureEdgeSmoothingOff();
	//sincSmoother->Update();

	vtkNew<vtkSTLWriter> stlWriter;
	stlWriter->SetFileName(fileName.c_str());
	//stlWriter->SetInputConnection(sincSmoother->GetOutputPort());
	stlWriter->SetInputConnection(norms->GetOutputPort());
	stlWriter->Write();
};

// -----------------------------------------------------------------
// Scaffold Face Builder

// constructor

