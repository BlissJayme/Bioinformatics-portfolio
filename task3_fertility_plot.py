# Python portfolio: Task 3 

# The code reads the fertility data from a CSV file, extracts the data from the selected countries and plots the average number of children per woman over time.

# Import the plotting library from matplotlib
import matplotlib.pyplot as plt

# The code reads the fertility data from the CSV file into a dictionary
# Use 'with open' to open the file
# The code reads the lines and saves to a list
with open("children-born-per-woman.csv") as f: 
    lines = f.readlines()

# The data is stored using a dictionary
# For each country, it contain a list of years and fertility rates
fertility_data = {}

# The code loops through every line in the list and skips the header line
# The functions .strip() removes the \n character and .split() splits the line into a list of items using the comma to separate the items
for line in lines[1:]:
    items = line.strip().split(",")

    

# Extract and convert the year and fertility rate
    country = items[0]
    year = int(items[2])
    fertility_rate = float(items[3])


# The code checks if the country has not been seen in the dictionary yet
    if country not in fertility_data:
        fertility_data[country] = {"Years": [], "Fertility_rates": []}

# The year and fertility rate are appended to the lists for this country
    fertility_data[country]["Years"].append(year)
    fertility_data[country]["Fertility_rates"].append(fertility_rate)



# Plotting the graph for the selected countries
# The names match the countries used in the CSV file
countries_to_plot = ['United Kingdom', 'Japan', 'China', 'Ecuador', 'Zimbabwe', 'India']

# Using matplotlib to plot the data
for country in countries_to_plot:
# The code plots the  country if it's present in the dictionary
    if country in fertility_data:
        data = fertility_data[country]

        #The years and fertility rates are taken directly from the dictionary for plotting
        years = data["Years"]
        rates = data["Fertility_rates"]

        # The fertility rate is plotted over time for the chosen country  
        plt.plot(years, rates, label=country)

# The code adds a legend to identify the countries
plt.legend()

# A title is added with the x - axis (year) and y - axis (fertility rate)
plt.title("Fertility rates over time in UK, Japan, China, Ecuador, Zimbabwe, and India")
plt.xlabel("Year")
plt.ylabel("Fertility rate (children per woman)")

# Limits are added for both the x and y axis (years 1950–2023 and fertility rate between 0 and 10)
plt.xlim(1950, 2023)
plt.ylim(0, 10)

# Display the graph
plt.show()

# Acknowledgement of AI use - I used ChatGPT to assist me with the location of the CSV file in my personal directory.