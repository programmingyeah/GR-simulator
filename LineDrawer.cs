using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class LineDrawer : MonoBehaviour
{

    [SerializeField]
    public Transform body;

    [SerializeField]
    public GameObject circlePrefab;

    [SerializeField]
    public GameObject simulation;

    private GameObject[] circleTrajectory;

    public LineRenderer lineRenderer;
    private Vector3[] positions;

    private int mode = 1;

    // Start is called before the first frame update
    void Start()
    {
        circleTrajectory = new GameObject[8000];
        positions = new Vector3[16000];

        if (lineRenderer == null)
        {
            lineRenderer = GetComponent<LineRenderer>();
        }

        // Set the positions for the LineRenderer
        lineRenderer.positionCount = positions.Length;
    }

    int currentPositionIndex = 0;
    // Update is called once per frame
    void Update()
    {
        if (mode == 0)
        {
            positions[currentPositionIndex] = body.localPosition;
            lineRenderer.positionCount = positions.Length;
            lineRenderer.SetPositions(positions);
            currentPositionIndex = (currentPositionIndex + 1) % 15999;
            UnityEngine.Debug.Log(positions[0]);
        }
        else if (mode == 1)
        {
            circleTrajectory[currentPositionIndex] = Instantiate(circlePrefab, body);
            circleTrajectory[currentPositionIndex].GetComponent<SpriteRenderer>().color = new Color(255,255,255);
            circleTrajectory[currentPositionIndex].transform.SetParent(simulation.transform);

            currentPositionIndex = (currentPositionIndex + 1) % 7999;
        }
    }
}
