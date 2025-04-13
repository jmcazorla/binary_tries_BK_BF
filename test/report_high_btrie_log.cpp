#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <iomanip>
#include <stdexcept>
#include <limits>

struct OutputData {
    size_t hsSize;
    double timeRatio;
    size_t intersectionSize;
    size_t clusterCounts[4];
    uint64_t minVal;
    uint64_t maxVal;
    double mean;
    double median;
    double stddev;
    double clusterLimit[3];
};

template <typename T>
double calculateAttributeMean(const std::vector<OutputData>& data, T OutputData::* memberPtr) {
    if (data.empty()) {
        return std::numeric_limits<double>::quiet_NaN(); // Return NaN for empty data
    }

    double sum = 0.0;
    for (const auto& item : data) {
        sum += static_cast<double>(item.*memberPtr);
    }

    double mean = sum / data.size();
    return mean;
}

// Function to group OutputData elements by timeRatio ranges using the mean
std::vector<std::vector<OutputData>> groupByTimeRatioUsingMean(const std::vector<OutputData>& allQueryData) {
    std::vector<std::vector<OutputData>> groups(4);
    if (allQueryData.empty()) {
        return groups;
    }

    double meanTime = calculateAttributeMean(allQueryData, &OutputData::timeRatio);

    // Define ranges based on the mean (you can adjust these factors)
    double range1 = meanTime / 2.0;
    double range2 = meanTime;
    double range3 = meanTime * 1.5; // Ejemplo: 1.5 veces la media

    for (const auto& data : allQueryData) {
        if (data.timeRatio <= range1) {
            groups[0].push_back(data);
        } else if (data.timeRatio > range1 && data.timeRatio <= range2) {
            groups[1].push_back(data);
        } else if (data.timeRatio > range2 && data.timeRatio <= range3) {
            groups[2].push_back(data);
        } else {
            groups[3].push_back(data);
        }
    }
    return groups;
}

std::pair<uint64_t, uint64_t> countUpDownAverage(const std::vector<OutputData>& data, double average) {
    uint64_t countUp = 0;
    uint64_t countDown = 0;
    for (const auto& item : data) {
        if (item.minVal > average) {
            countUp += item.clusterCounts[0];
        }
        if (item.clusterLimit[0] > average) {
            countUp += item.clusterCounts[1];
        }
        if (item.clusterLimit[1] > average) {
            countUp += item.clusterCounts[2];
        }
        if (item.clusterLimit[2] > average) {
            countUp += item.clusterCounts[3];
        }

        if (item.clusterLimit[0] < average) {
            countDown += item.clusterCounts[0];
        }
        if (item.clusterLimit[1] < average) {
            countDown += item.clusterCounts[1];
        }
        if (item.clusterLimit[2] < average) {
            countDown += item.clusterCounts[2];
        }
        if (item.maxVal < average) {
            countDown += item.clusterCounts[3];
        }
    }
    return {countUp/data.size(), countDown/data.size()};
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <binary_file_name>" << std::endl;
        return 1;
    }

    const std::string binaryFileName = argv[1];
    std::ifstream inputFile(binaryFileName, std::ios::binary);
    std::vector<OutputData> allQueryData;

    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open binary file '" << binaryFileName << "' for reading." << std::endl;
        return 1;
    }

    OutputData currentData;
    while (inputFile.read(reinterpret_cast<char*>(&currentData), sizeof(OutputData))) {
        allQueryData.push_back(currentData);
    }
    inputFile.close();

    if (allQueryData.empty()) {
        std::cout << "No OutputData records found in the file." << std::endl;
        return 0;
    }

    std::cout << "\nReport of Mean Values for OutputData Attributes:\n";
    std::cout << "---------------------------------------------------\n";
    double mean = 0.0;
    double meanTime = 0.0;
    double average = 0.0;
    mean = calculateAttributeMean(allQueryData, &OutputData::hsSize);
    std::cout << "Average query size : " << std::fixed << std::setprecision(2) << mean << std::endl;
    meanTime = calculateAttributeMean(allQueryData, &OutputData::timeRatio);
    std::cout << "Average time: " << std::fixed << std::setprecision(2) << meanTime << std::endl;
    mean = calculateAttributeMean(allQueryData, &OutputData::intersectionSize);
    std::cout << "Average number of subsets: " << std::fixed << std::setprecision(2) << mean << std::endl;
    mean = calculateAttributeMean(allQueryData, &OutputData::minVal);
    std::cout << "Average minimum subset sizes: " << std::fixed << std::setprecision(2) << mean << std::endl;
    mean = calculateAttributeMean(allQueryData, &OutputData::maxVal);
    std::cout << "Average maximum subset sizes: " << std::fixed << std::setprecision(2) << mean << std::endl;
    average = calculateAttributeMean(allQueryData, &OutputData::mean);
    std::cout << "Average subset size average: " << std::fixed << std::setprecision(2) << average << std::endl;
    mean = calculateAttributeMean(allQueryData, &OutputData::median);
    std::cout << "Average of the medians of the subset sizes: " << std::fixed << std::setprecision(2) << mean << std::endl;
    mean = calculateAttributeMean(allQueryData, &OutputData::stddev);
    std::cout << "Average standard deviation of subset size: " << std::fixed << std::setprecision(2) << mean << std::endl;

    std::cout << "---------------------------------------------------\n";

    std::vector<std::vector<OutputData>> groupByTime = groupByTimeRatioUsingMean(allQueryData);
    std::cout << "Grouped by time ratio using mean:\n";

    double range[3] = {meanTime / 2.0, meanTime, meanTime * 1.5}; // Ejemplo: 1.5 veces la media

    for (size_t i = 0; i < groupByTime.size(); ++i) {
        std::cout << "Group " << i + 1 << ": " << groupByTime[i].size() << " records\n";
        std::cout << "Time ratio range: ";
        if (i == 0) {
            std::cout << " <= " << range[0];
        } else if (i == 1) {
            std::cout << " > " << range[0] << " and <= " << range[1];
        } else if (i == 2) {
            std::cout << " > " << range[1] << " and <= " << range[2];
        } else {
            std::cout << " > " << range[2];
        }
        std::cout << "\n";
        mean = 0.0;
        mean = calculateAttributeMean(groupByTime[i], &OutputData::hsSize);
        std::cout << "Average query size : " << std::fixed << std::setprecision(2) << mean << std::endl;
        mean = calculateAttributeMean(groupByTime[i], &OutputData::timeRatio);
        std::cout << "Average time: " << std::fixed << std::setprecision(2) << mean << std::endl;
        mean = calculateAttributeMean(groupByTime[i], &OutputData::intersectionSize);
        std::cout << "Average number of subsets: " << std::fixed << std::setprecision(2) << mean << std::endl;
        mean = calculateAttributeMean(groupByTime[i], &OutputData::minVal);
        std::cout << "Average minimum subset sizes: " << std::fixed << std::setprecision(2) << mean << std::endl;
        mean = calculateAttributeMean(groupByTime[i], &OutputData::maxVal);
        std::cout << "Average maximum subset sizes : " << std::fixed << std::setprecision(2) << mean << std::endl;
        mean = calculateAttributeMean(groupByTime[i], &OutputData::mean);
        std::cout << "Average subset size average: " << std::fixed << std::setprecision(2) << mean << std::endl;
        mean = calculateAttributeMean(groupByTime[i], &OutputData::median);
        std::cout << "Average of the medians of the subset sizes: " << std::fixed << std::setprecision(2) << mean << std::endl;
        mean = calculateAttributeMean(groupByTime[i], &OutputData::stddev);
        std::cout << "Average standard deviation of subset size: " << std::fixed << std::setprecision(2) << mean << std::endl;

        std::cout << "---------------------------------------------------\n";
    }


    return 0;
}