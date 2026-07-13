# UAV Multi-Sensor Fusion and Kalman Filter Estimation 🚁

This project implements a 9-axis orientation (Roll, Pitch, Yaw) and altitude estimation system for Unmanned Aerial Vehicles (UAVs) using a multi-sensor fusion architecture.

## 📌 Project Overview
The system collects raw data from MPU6050, QMC5883P, and BMP280 sensors using an **STM32F072B** microcontroller. The data is transmitted via UART and processed offline in **MATLAB** using a custom discrete-time Kalman Filter to eliminate high-frequency noise, gyroscope drift, and magnetic interference.

## 📂 Project Structure
* `/MATLAB_and_Sensor_Data`: Contains the raw sensor dataset (`sensor_data.csv`) and the MATLAB script (`Data_Reading.m`) for Kalman Filter processing and data analysis.
* `/Reports`: Detailed technical project reports.
  * [📄 Read the English Report here](Reports/STM32_MATLAB_Kalman_Filter_English_Report.pdf)
  * [📄 Türkçe Raporu buradan okuyabilirsiniz](Reports/STM32_MATLAB_Kalman_Filter_Türkçe_Rapor.pdf)
* `STM32.zip`: Compressed archive containing the embedded C codes, libraries, and STM32CubeIDE project files.

## 🛠️ Hardware & Software Used
* **Hardware:** STM32F072B Discovery, MPU6050, QMC5883P, BMP280
* **Software/Tools:** STM32CubeIDE (C), MATLAB, Tera Term

## 👤 Author
* **Süleyman Açıkal**
