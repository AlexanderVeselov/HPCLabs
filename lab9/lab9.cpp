/*

2. Первая задача про Китайский банк. На маленькой улице Чжуань-
Сю в городе Гонконг живут двести тысяч китайцев и находятся три банка.
Каждый из этих банков принимает деньги от вкладчиков в трех валютах –
китайских юанях, американских долларах и английских фунтах стерлингов.
При этом если вкладчик хочет взять деньги в одном банке на улице Чжуань-
Сю и положить в другой, то ему в первом банке выдается только расписка,
которую он и относит во второй банк. В пятницу вечером банки подсчиты-
вают, сколько денег и в какой валюте они должны соседям и отправляют ин-
кассаторов отнести эти деньги. Написать программу, моделирующую обмен
деньгами в пятницу вечером на улице Чжуань-Сю, используя метод передачи
информации «точка-точка».
3. Вторая задача про Китайский банк. Решить задачу о китайских
банках с дополнительным условием, что все пятничные расчеты банки про-
водят через Большой Банк, находящийся на улице Чжуань-Го.

*/

#include <mpi.h>
#include <iostream>
#include <unordered_map>
#include <array>
#include <algorithm>
#include <numeric>

enum Currency
{
    kDollar = 0,
    kYuan,
    kPound,

    kNumCurrencies
};

static std::unordered_map<Currency, float> g_exchange_rates =
{
    { Currency::kDollar, 1.0f },
    { Currency::kYuan,   0.14f },
    { Currency::kPound,  1.28f }
};

class Money
{
public:
    Money() {}
    Money(float amount, Currency currency)
        : amount_(amount), currency_(currency)
    {}

    float GetAmount() const
    {
        return amount_;
    }

    Currency GetCurrency() const
    {
        return currency_;
    }

    Money ToCurrency(Currency new_currency) const
    {
        return Money(amount_ * g_exchange_rates[currency_] / g_exchange_rates[new_currency], new_currency);
    }

    Money& operator+=(Money const& other)
    {
        amount_ += (other.currency_ == currency_) ? other.amount_ : other.ToCurrency(currency_).amount_;
        return *this;
    }

    Money& operator-=(Money const& other)
    {
        amount_ -= (other.currency_ == currency_) ? other.amount_ : other.ToCurrency(currency_).amount_;
        return *this;
    }

private:
    float amount_;
    Currency currency_;
};

struct TransferTicket
{
    Money money;
    // int - to avoid additional casts
    int dst_bank_id;
    std::size_t src_account_id;
    std::size_t dst_account_id;

};

class CentralBank;

class Bank
{
public:
    Bank();
    void GenerateRandomTickets();
    void PrintBankAccounts() const;
    void PrintTickets() const;
    int GetId() const { return bank_id_; }
    void SendMoney();
    const int kCentralBankId = 0;

private:
    int num_banks_;
    int bank_id_;
    std::vector<Money> bank_accounts_;
    std::vector<TransferTicket> tickets_;

};

char GetCurrencySymbol(Currency currency)
{
    static const std::unordered_map<Currency, char> currency_symbols =
    {
        { Currency::kDollar, '$' },
        { Currency::kYuan,   'Y' },
        { Currency::kPound,  'E' }
    };

    return currency_symbols.find(currency)->second;
}

std::ostream & operator<< (std::ostream & os, Money money)
{
    os << money.GetAmount() << " " << GetCurrencySymbol(money.GetCurrency());
    return os;
}

Bank::Bank()
{
    // Get current bank id
    MPI_Comm_rank(MPI_COMM_WORLD, &bank_id_);
    std::srand(bank_id_);

    MPI_Comm_size(MPI_COMM_WORLD, &num_banks_);

    std::size_t num_accounts = 3;
    // Generate random accounts
    bank_accounts_.resize(num_accounts);
    std::generate(bank_accounts_.begin(), bank_accounts_.end(), []()
    {
        return Money((float)std::rand() / RAND_MAX * 100.0f, Currency::kDollar);//(Currency)(std::rand() % Currency::kNumCurrencies));
    });

}

void Bank::GenerateRandomTickets()
{
    if (bank_id_ == kCentralBankId)
    {
        return;
    }

    // Generate random tickets
    std::size_t num_tickets = std::rand() % 3;
    tickets_.resize(num_tickets);
    std::generate(tickets_.begin(), tickets_.end(), [this]()
    {
        int dst_bank_id;
        do
        {
            dst_bank_id = std::rand() % num_banks_;
        } while (dst_bank_id == bank_id_ || dst_bank_id == kCentralBankId);

        TransferTicket ticket
        {
            Money((float)std::rand() / RAND_MAX * 10.0f, Currency::kDollar),//(Currency)(std::rand() % Currency::kNumCurrencies)),
            dst_bank_id,
            // src id
            std::rand() % bank_accounts_.size(),
            // dst id
            std::rand() % 3
        };

        return ticket;

    });
}

void Bank::PrintBankAccounts() const
{
    for (std::size_t i = 0; i < bank_accounts_.size(); ++i)
    {
        std::cout << "id: " << i << " balance: " << bank_accounts_[i] << std::endl;
    }
}

void Bank::PrintTickets() const
{
    for (std::size_t i = 0; i < tickets_.size(); ++i)
    {
        std::cout << "Transfer " << tickets_[i].money << " to bank " << tickets_[i].dst_bank_id << " from id "
            << tickets_[i].src_account_id << " to id " << tickets_[i].dst_account_id << std::endl;
    }

}

//template <class T>
//void GatherData(std::vector<T> const& src_data, std::vector<T> & dst_data, int master_id)
//{
//    int node_id;
//    MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
//
//    std::array<int, 4> recv_counts;
//    int send_counts = src_data.size();
//    MPI_Gather(&send_counts, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, master_id, MPI_COMM_WORLD);
//
//    std::size_t total_counts;
//    if (node_id == master_id)
//    {
//        total_counts = std::accumulate(recv_counts.begin() + 1, recv_counts.end(), 0);
//    }
//}

template <class T>
void SendPointPoint(std::vector<T> const& src_data, std::vector<T> & dst_data, int src_node_id, int dst_node_id)
{
    int node_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &node_id);

    int count = src_data.size();

    // Send count
    if (node_id == src_node_id)
    {
        MPI_Send(&count, 1, MPI_INT, dst_node_id, 0, MPI_COMM_WORLD);
    }

    // Recv count
    if (node_id == dst_node_id)
    {
        MPI_Status status;
        MPI_Recv(&count, 1, MPI_INT, src_node_id, 0, MPI_COMM_WORLD, &status);

        dst_data.clear();
        //std::cout << "Node " << dst_node_id << " - resize dst_data to " << count << std::endl;
        dst_data.resize(count);

    }

    // Send data
    if (node_id == src_node_id)
    {
        MPI_Send(src_data.data(), count * sizeof(T), MPI_CHAR, dst_node_id, 0, MPI_COMM_WORLD);
    }

    // Recv data
    if (node_id == dst_node_id)
    {
        MPI_Status status;
        MPI_Recv(dst_data.data(), count * sizeof(T), MPI_CHAR, src_node_id, 0, MPI_COMM_WORLD, &status);
    }

}

void Bank::SendMoney()
{
    // Write off local money
    for (std::size_t i = 0; i < tickets_.size(); ++i)
    {
        TransferTicket const& ticket = tickets_[i];
        bank_accounts_[ticket.src_account_id] -= ticket.money;

    }

    // Calculate total tickets count
    std::array<int, 4> recv_ticket_counts;
    int send_ticket_counts = tickets_.size();
    MPI_Gather(&send_ticket_counts, 1, MPI_INT, recv_ticket_counts.data(), 1, MPI_INT, kCentralBankId, MPI_COMM_WORLD);

    std::size_t total_tickets;
    if (bank_id_ == kCentralBankId)
    {
        total_tickets = std::accumulate(recv_ticket_counts.begin() + 1, recv_ticket_counts.end(), 0);
    }

    // Received data
    std::vector<TransferTicket> central_tickets;

    if (bank_id_ == kCentralBankId)
    {
        central_tickets.resize(total_tickets);
    }

    std::array<int, 4> recv_byte_counts;
    for (std::size_t i = 0; i < 4; ++i)
    {
        recv_byte_counts[i] = recv_ticket_counts[i] * sizeof(TransferTicket);
    }

    std::array<int, 4> recv_byte_offsets;
    recv_byte_offsets[0] = 0;
    for (std::size_t i = 1; i < 4; ++i)
    {
        recv_byte_offsets[i] = recv_byte_offsets[i - 1] + recv_byte_counts[i - 1];
    }

    MPI_Gatherv(tickets_.data(), tickets_.size() * sizeof(TransferTicket), MPI_CHAR, central_tickets.data(), recv_byte_counts.data(), recv_byte_offsets.data(), MPI_CHAR, kCentralBankId, MPI_COMM_WORLD);

    std::vector<TransferTicket> recv_tickets;

    for (std::size_t dst_bank = 1; dst_bank < num_banks_; ++dst_bank)
    {
        if (bank_id_ == kCentralBankId)
        {
            std::copy_if(central_tickets.begin(), central_tickets.end(), std::back_inserter(recv_tickets), [dst_bank](TransferTicket const& ticket)
            {
                return ticket.dst_bank_id == dst_bank;
            });
        }

        SendPointPoint(recv_tickets, recv_tickets, kCentralBankId, dst_bank);

        if (bank_id_ == kCentralBankId)
        {
            recv_tickets.clear();
        }
    }

    // Top up accounts
    if (bank_id_ != kCentralBankId)
    {
        for (std::size_t i = 0; i < recv_tickets.size(); ++i)
        {
            TransferTicket const& ticket = recv_tickets[i];
            bank_accounts_[ticket.dst_account_id] += ticket.money;

        }
    }

}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    Bank bank;
    bank.GenerateRandomTickets();

    for (int i = 1; i <= 3; ++i)
    {
        if (bank.GetId() == i)
        {
            std::cout << "Bank " << i << " accounts: " << std::endl;
            bank.PrintBankAccounts();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    for (int i = 1; i <= 3; ++i)
    {
        if (bank.GetId() == i)
        {
            std::cout << "Bank " << i << " tickets: " << std::endl;
            bank.PrintTickets();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    bank.SendMoney();


    for (int i = 1; i <= 3; ++i)
    {
        if (bank.GetId() == i)
        {
            std::cout << "Bank " << i << " accounts: " << std::endl;
            bank.PrintBankAccounts();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
