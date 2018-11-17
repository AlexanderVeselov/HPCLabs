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
#include <deque>
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
    std::deque<TransferTicket> tickets_;

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

    std::size_t num_accounts = 20;
    // Generate random accounts
    bank_accounts_.resize(num_accounts);
    std::generate(bank_accounts_.begin(), bank_accounts_.end(), []()
    {
        return Money((float)std::rand() / RAND_MAX * 100.0f, (Currency)(std::rand() % Currency::kNumCurrencies));
    });

}

void Bank::GenerateRandomTickets()
{
    if (bank_id_ == kCentralBankId)
    {
        return;
    }

    // Generate random tickets
    std::size_t num_tickets = std::rand() % 20;
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
            Money((float)std::rand() / RAND_MAX * 10.0f, (Currency)(std::rand() % Currency::kNumCurrencies)),
            dst_bank_id,
            // src id
            std::rand() % bank_accounts_.size(),
            // dst id
            std::rand() % bank_accounts_.size()
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

void Bank::SendMoney()
{
    std::vector<float> send_amounts;
    std::vector<Currency> send_currencies;
    std::vector<int> send_dst_bank_ids;
    std::vector<std::size_t> send_dst_account_ids;

    if (bank_id_ != kCentralBankId)
    {
        send_amounts.resize(tickets_.size());
        send_currencies.resize(tickets_.size());
        send_dst_bank_ids.resize(tickets_.size());
        send_dst_account_ids.resize(tickets_.size());

        // Serialize tickets
        for (std::size_t i = 0; i < tickets_.size(); ++i)
        {
            TransferTicket const& ticket = tickets_[i];
            // Write off money
            bank_accounts_[ticket.src_account_id] -= ticket.money;

            // Serialize ticket
            send_amounts[i] = ticket.money.GetAmount();
            send_currencies[i] = ticket.money.GetCurrency();
            send_dst_bank_ids[i] = ticket.dst_bank_id;
            send_dst_account_ids[i] = ticket.dst_account_id;
        }
    }

    // Calculate total tickets count
    std::array<int, 4> recv_counts;
    std::size_t send_tickets = tickets_.size();
    MPI_Gather(&send_tickets, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, kCentralBankId, MPI_COMM_WORLD);

    std::size_t total_tickets;
    if (bank_id_ == kCentralBankId)
    {
        total_tickets = std::accumulate(recv_counts.begin() + 1, recv_counts.end(), 0);
        std::cout << "Recv counts " << total_tickets << std::endl;
        for (int i = 0; i < recv_counts.size(); ++i)
        {
            std::cout << recv_counts[i] << std::endl;
        }
    }

    std::vector<float> recv_amounts;
    std::vector<Currency> recv_currencies;
    std::vector<int> recv_dst_bank_ids;
    std::vector<std::size_t> recv_dst_account_ids;

    if (bank_id_ == kCentralBankId)
    {
        recv_amounts.resize(total_tickets);
        recv_currencies.resize(total_tickets);
        recv_dst_bank_ids.resize(total_tickets);
        recv_dst_account_ids.resize(total_tickets);
    }

    std::array<int, 4> recv_offsets;
    recv_offsets[0] = 0;
    for (std::size_t i = 1; i < 4; ++i)
    {
        recv_offsets[i] = recv_offsets[i - 1] + recv_counts[i - 1];
    }

    if (bank_id_ == kCentralBankId)
    {
        std::cout << "Recv offsets " << recv_offsets.size() << std::endl;
        for (int i = 0; i < recv_offsets.size(); ++i)
        {
            std::cout << recv_offsets[i] << std::endl;
        }
    }

    MPI_Gatherv(send_amounts.data(), send_amounts.size(), MPI_FLOAT, recv_amounts.data(), recv_counts.data(), recv_counts.data(), MPI_FLOAT, kCentralBankId, MPI_COMM_WORLD);

    //if (bank_id_ == kCentralBankId)
    //{
    //    std::cout << "Send amounts " << recv_amounts.size() << std::endl;
    //    for (int i = 0; i < recv_amounts.size(); ++i)
    //    {
    //        std::cout << recv_amounts[i] << std::endl;
    //    }
    //}


}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    Bank bank;
    bank.GenerateRandomTickets();
    bank.SendMoney();

    MPI_Finalize();
    return 0;
}
